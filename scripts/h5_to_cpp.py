#!/usr/bin/env python3

import argparse
import builtins
import contextlib
import itertools
import os.path
import sys

import h5py
import numpy as np

# -------------------------------------------------------

numpy_to_hdf5 = {
    "i": "STD_I",
    "l": "STD_I",
    "B": "STD_U",
    "L": "STD_I",
    "f": "IEEE_F",
    "d": "IEEE_F",
}

# -------------------------------------------------------

hdr_template = \
"""
//  WARNING:
//    This file was autogenerated by {filename}.
//    Do not edit by hand!
//    
//    The invocation that generated this file was:
//
//       {cmdline}
//


#ifndef {guard_var}
#define {guard_var}

#include <array>

#include "H5Cpp.h"
#include "readH5/BufferView.h"

namespace {namespace}
{{

  // This generic template will be overloaded
  // for every specific type that we create below
  template <typename T>
  H5::CompType BuildCompType();

{members}

}}

#endif // {guard_var}
"""

# -----------

impl_template = \
"""
//  WARNING:
//    This file was autogenerated by {filename}.
//    Do not edit by hand!
//    
//    The invocation that generated this file was:
//
//       {cmdline}
//

#include "{hdr}"
#include "H5Cpp.h"

namespace {namespace}
{{

{members}

}}
"""

# -----------

stdlib_include_template = "#include <{header}>"

# -----------

class_template = \
"""
struct {name}
{{
{members}
}};
"""

# -----------

simple_member_template = "{typ} {name};"

# -----------

region_ref_fn_template = """
template <typename T>
const hdset_reg_ref_t& GetRef() const
{{
{members}
}}
"""

# -----------

region_ref_member_template = "{if_expr}(std::is_same_v<T, {typename}>) return {fieldname};"

# -----------

handle_members_note = \
"""
// note: the following 'handle' objects
// are used internally by HDF5 to keep track
// of the memory for variable-length buffers.
// please use the SyncVectors() method
// after loading data into the object
// to fill the corresponding BufferView<>s above,
// and then use those for access to the data.
"""

# -----------

sync_method_template = \
"""
void {klass}::SyncVectors()
{{
{members}
}}
"""

# -----------

sync_method_member_template = "{var}.reset(&{handle});"

# -----------

fwd_declare_template = "struct {typ};"

# -----------


enum_cpp_template = \
"""
enum class {name} : {size_type}
{{
{members}
}};
"""

# -----------

enum_cpp_member_template = "{name} = {val},"

# -----------

compound_type_decl_template = \
"""
template <>
H5::CompType BuildCompType<{klass}>();
"""

# -----------

compound_type_impl_template = \
"""
template <>
H5::CompType BuildCompType<{klass}>()
{{
  H5::CompType ctype(sizeof({klass}));

{members}

  return ctype;
}}
"""

# -----------

compound_type_member_template = 'ctype.insertMember("{h5_name}", HOFFSET({klass}, {cpp_name}), {h5type});'

# -----------

compound_type_enum_member_template = \
"""
H5::EnumType {name}_enumtype({h5type});
{enum_type} {name}_enum_val;
{members}
ctype.insertMember("{name}", HOFFSET({klass}, {name}), {name}_enumtype);
"""

# -----------

compound_type_enum_entry_template = '{name}_enum_val = {val}; {name}_enumtype.insert("{h5name}", &{name}_enum_val);'

# -----------

compound_type_string_template = \
"""
H5::StrType {h5_name}_strType({strtype}, {len});
{h5_name}_strType.setCset({charset});
ctype.insertMember("{h5_name}", HOFFSET({klass}, {cpp_name}), {h5_name}_strType);
"""

# -------------------------------------------------------

def dataset_to_name(dataset):
    assert hasattr(dataset, "name")
    return dataset.name.lstrip("/")

# -------------------------------------------------------


class Serializable:
    """ Takes a string template and (recursively, if necessary) fills it in with strings from its members """
    def __init__(self, template, template_args, member_list=(), member_join="\n", member_indent="  ", base_indent=""):
        self.template = template
        self.template_args = template_args
        self.members = member_list
        self.join_str = member_join
        self.member_indent = member_indent
        self.base_indent = base_indent

        assert all(hasattr(member, "emit") for member in member_list), "Incompatible member type in list: " + str(
            member_list)

    def emit(self, indent=""):
        member_str = self.join_str.join(
            self.base_indent + member.emit(indent=self.member_indent) for member in self.members)

        args = self.template_args.copy()
        args["members"] = member_str

        out = self.template.format(**args)
        return "\n".join([indent + self.base_indent + line for line in out.split("\n")])

# -------------------------------------------------------

class TypeSerializer:
    def __init__(self, class_name_map):
        self.class_name_map = class_name_map

        self.discovered_enums = {}
        self.enum_vals_h5 = {}
        self.cpp_types = {}
        self.cpp_types_impl = {}
        self.comptype_builders_decl = {}
        self.comptype_builders_impl = {}
        self.fwd_declares = {}

        self.cpp_headers = []

        # whenever we need to regenerate the Serializable
        self._dirty = False
        self._serializables = {}

    def type_string(self, typ, fieldname=None, which="cpp"):
        assert which in ("cpp", "h5"), "Unrecognized kind of type string: " + which

        # sys.byteorder uses these instead of '<' and '>' like numpy or 'LE' and 'BE' like hdf5
        sysorder = {
            "little": "LE",
            "big": "BE",
        }

        cpp_name = ""
        h5_name = ""

        if not isinstance(typ, builtins.type) and h5py.check_dtype(ref=typ):
            cpp_name += "hdset_reg_ref_t"
            h5_name = "H5::PredType::STD_REF_DSETREG"
        # horrible hack, but bools are always stored as ints, so we have no other way of knowing
        elif fieldname and fieldname.startswith("is_"):
            cpp_name += "bool"
            h5_name = self.type_string(typ, fieldname=None, which="h5")
        elif h5py.check_enum_dtype(typ):
            typenames = fieldname.split("_")
            typename = "".join(t.capitalize() for t in typenames)

            # we really shouldn't have the side effect of storing the enums
            # performed inside a function that's doing something else,
            # but it's the simplest way to do it
            self.discovered_enums[typename] = Serializable(template=enum_cpp_template, template_args={"name": typename,
                                                                                                      "size_type" : self.type_string(np.dtype(str(typ)))},
                                                           member_list=[Serializable(template=enum_cpp_member_template,
                                                                                     template_args=dict(
                                                                                         name="k{name}".format(
                                                                                             name="".join(k.split())),
                                                                                         val=v))
                                                                        for k, v in h5py.check_enum_dtype(typ).items()])
            self.enum_vals_h5[fieldname] = {k: Serializable(template=compound_type_enum_entry_template,
                                                            template_args=dict(name=fieldname, val=v, h5name=k))
                                            for k, v in h5py.check_enum_dtype(typ).items()}

            h5_name = self.type_string(np.dtype(str(typ)), fieldname, which="h5")
            # print("discovered enum:", typename, self.discovered_enums[typename])
            cpp_name += typename
        elif typ == str or h5py.check_string_dtype(typ):
            cpp_name += "char *"
            h5_name = "H5::PredType::H5T_C_S1"
        elif h5py.check_vlen_dtype(typ):
            cpp_name += "BufferView<{base_type}>".format(base_type=self.type_string(h5py.check_vlen_dtype(typ), fieldname))
            h5_name = "H5::VarLenType({typ})".format(typ=self.type_string(h5py.check_vlen_dtype(typ), fieldname, which="h5"))
        elif typ.ndim > 0:
            assert typ.ndim == 1, "Don't know how to handle multi-dimensional array types"
            cpp_name += "std::array<{type}, {len}>".format(type=self.type_string(typ.subdtype[0], fieldname),
                                                           len=typ.shape[0])
            h5_name = "H5::ArrayType({typ}, 1, &std::array<hsize_t, 1>{{{count}}}[0])".format(
                typ=self.type_string(typ.subdtype[0], fieldname, which="h5"), count=typ.shape[0])
        elif typ.char in numpy_to_hdf5:
            h5_name = "H5::PredType::"
            h5_name += numpy_to_hdf5[typ.char]
            h5_name += str(typ.itemsize * 8)  # numpy reports bytes, HDF5 uses bits
            h5_name += "BE" if typ.byteorder == ">" else "LE" if typ.byteorder == "<" else sysorder[sys.byteorder]

            if typ.name.startswith("int") or typ.name.startswith("uint"):
                cpp_name = typ.name + "_t"  # e.g.: int64_t
            elif typ.name.startswith("float"):
                if typ.name == "float32":
                    cpp_name = "float"
                elif typ.name == "float64":
                    cpp_name = "double"
                elif typ.name == "float128":
                    cpp_name = "long double"
            assert len(cpp_name) > 0, "Couldn't understand type name: '{0}'".format(typ.name)

        else:
            print("typ.char =", typ.char)
            raise TypeError("Don't know how to handle type: " + str(typ))

        h5_name = h5_name.replace("H5T_", "")

        return cpp_name if which == "cpp" else h5_name

    def add_dataset(self, dataset):
        ds_name = dataset_to_name(dataset)
        assert ds_name in self.class_name_map, "Unrecognized dataset name: '{0}'".format(ds_name)
        class_name = self.class_name_map[ds_name]

        self._dirty = True

        cpp_members = []
        h5_members = []
        handles = []
        region_refs = {}
        print("Examining dataset:", dataset.name)
        for fieldname in dataset.dtype.names:
            typ = dataset.dtype[fieldname]
            cpp_members.append(Serializable(template=simple_member_template,
                                            template_args=dict(name=fieldname,
                                                               typ=self.type_string(typ,fieldname=fieldname))))

            # for variable-length items, we have to maintain a special "handle"
            # in addition to the vector generated for cpp_members
            # (unfortunately check_vlen_dtype() will return True for string types...)
            if not h5py.check_string_dtype(typ) and h5py.check_vlen_dtype(typ):
                handles.append(fieldname)

            # for region references, we need to connect the dtype to the field name
            if h5py.check_dtype(ref=typ):
                ref_name = dataset_to_name(dataset.file[fieldname])
                if ref_name in self.class_name_map:
                    region_refs[self.class_name_map[ref_name]] = fieldname

            # enums need to be handled specially because they have a different template
            if h5py.check_enum_dtype(typ):
                h5_members.append(Serializable(template=compound_type_enum_member_template,
                                               template_args=dict(name=fieldname,
                                                                  klass=class_name,
                                                                  enum_type=self.type_string(np.dtype(str(typ)), fieldname=fieldname),
                                                                  h5type=self.type_string(typ, fieldname=fieldname, which="h5"),
                                                                  ),
                                               member_list=self.enum_vals_h5[fieldname].values()))
            # strings too
            elif h5py.check_string_dtype(typ):
                info = h5py.check_string_dtype(typ)
                h5_members.append(Serializable(template=compound_type_string_template,
                                               template_args=dict(h5_name=fieldname,
                                                                  cpp_name=fieldname,
                                                                  klass=class_name,
                                                                  len="H5T_VARIABLE" if info.length is None else info.length,
                                                                  charset="H5T_CSET_UTF8" if info.encoding == "utf-8" else "H5T_CSET_ASCII",
                                                                  strtype=self.type_string(typ, fieldname=fieldname, which="h5"),)
                                               )
                                  )
            else:
                h5_members.append(Serializable(template=compound_type_member_template,
                                               template_args=dict(h5_name=fieldname,
                                                                  cpp_name=fieldname + ("_handle" if h5py.check_vlen_dtype(typ) else ""),
                                                                  klass=class_name,
                                                                  h5type=self.type_string(typ, fieldname=fieldname, which="h5"))))
        sync_members = []
        if len(handles) > 0:
            sync_members = [Serializable(template=sync_method_member_template,
                                         template_args=dict(var=handle.rsplit("_handle", 1)[0],
                                                            typ=self.type_string(h5py.check_vlen_dtype(dataset.dtype[handle.rsplit("_handle", 1)[0]])),
                                                            handle=handle + "_handle")) for handle in handles]

        cpp_members.append(Serializable("\nvoid SyncVectors();", template_args={}))

        if len(handles) > 0:
            cpp_members.append(Serializable(template=handle_members_note, template_args={}))
            for handle in handles:
                cpp_members.append(Serializable(template=simple_member_template,
                                                template_args=dict(name=handle + "_handle", typ="hvl_t")))

        if len(region_refs) > 0:
            if_stmts = [Serializable(template=region_ref_member_template,
                                     template_args=dict(if_expr="if constexpr" if i==0 else "else if",
                                                        typename=typename,
                                                        fieldname=region_fieldname))
                        for i, (typename, region_fieldname) in enumerate(region_refs.items())]
            cpp_members.append(Serializable(template=region_ref_fn_template,
                                            template_args={},
                                            member_list=if_stmts))

        self.cpp_types[class_name] = Serializable(template=class_template, template_args=dict(name=class_name),
                                                  member_list=cpp_members)
        self.fwd_declares[class_name] = Serializable(template=fwd_declare_template, template_args=dict(typ=class_name))

        self.cpp_types_impl[class_name] = Serializable(template=sync_method_template,
                                                       template_args=dict(klass=class_name),
                                                       member_list=sync_members)
        self.comptype_builders_decl[class_name] = Serializable(template=compound_type_decl_template,
                                                               template_args=dict(klass=class_name))
        self.comptype_builders_impl[class_name] = Serializable(template=compound_type_impl_template,
                                                               template_args=dict(klass=class_name),
                                                               member_list=h5_members)

    def emit(self, namespace, header_filename):
        if len(self._serializables) == 0 or self._dirty:
            self._serializables[".h"] = Serializable(template=hdr_template,
                                                     template_args=dict(filename=os.path.basename(__file__),
                                                                        cmdline=" ".join(sys.argv),
                                                                        namespace=namespace,
                                                                        guard_var=namespace.replace("::", "_").upper() + "_" + header_filename.replace(".", "_").upper()),
                                                     member_list=tuple(itertools.chain(self.discovered_enums.values(),
                                                                                       self.fwd_declares.values(),
                                                                                       self.cpp_types.values(),
                                                                                       self.comptype_builders_decl.values())),
                                                     member_indent="  ")
            self._serializables[".cxx"] = Serializable(template=impl_template,
                                                       template_args=dict(filename=os.path.basename(__file__),
                                                                          cmdline=" ".join(sys.argv),
                                                                          namespace=namespace,
                                                                          hdr=header_filename),
                                                       member_list=tuple(itertools.chain(self.cpp_types_impl.values(),
                                                                                         self.comptype_builders_impl.values())),
                                                       member_indent="  ")

            self._dirty = False

        return {k : s.emit() for k, s in self._serializables.items()}

# -------------------------------------------------------
# -------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="Utility for autogenerating the necessary C++ code to read structured datatypes in an .h5 file.")

    parser.add_argument("-f", "--hdf5",
                        action="store",
                        required=True,
                        help="name of .h5 file to find dataset(s) in")

    parser.add_argument("-o", "--outfile_stub",
                        action="store",
                        required=True,
                        help="Base name of the .h and .cxx files that will be generated.  ('.h' and '.cxx' extensions will be added automatically.")

    parser.add_argument("-d", "--dataset",
                        action="append",
                        default=[],
                        help="Dataset(s) to examine (specify as many as needed).  If unspecified, will use all datasets in file")

    parser.add_argument("-cn", "--class_name",
                        action="append",
                        default=[],
                        help="Output class names for the datasets specified with --dataset.  " +
                             "Must have the same number of entries as --dataset.  " +
                             "If not specified, dataset names will be used (capitalized).")

    parser.add_argument("-ns", "--namespace",
                        default="test",
                        help="Namespace that generated C++ declarations & implementations should be nested under.")

    args = parser.parse_args()

    if len(args.dataset) != len(args.class_name):
        print("Datasets and class names do not match!", file=sys.stderr)
        sys.exit(1)

    datasets = args.dataset
    class_names = args.class_name

    f = h5py.File(args.hdf5, 'r')
    if len(args.dataset) == 0:
        datasets = [f[k] for k in f.keys() if hasattr(f[k], "dtype")]
        class_names = [ds.name.capitalize() for ds in datasets]
        assert len(datasets) > 0, "No datasets found in file: " + args.filename
#        print("found datasets:", datasets)
    else:
        datasets = []
        for ds in args.dataset:
            if ds not in f or not hasattr(f[ds], "dtype"):
                print("Dataset '{0}' not found in file: {1}".format(ds, args.filename), file=sys.stderr)
                sys.exit(1)
            datasets.append(f[ds])

    with contextlib.ExitStack() as stack:
        # this opens all of the files with a context manager
        # (so they all get closed properly, even if there's an exception)
        files = {suffix: stack.enter_context(open(args.outfile_stub + suffix, "w")) for suffix in (".h", ".cxx")}

        # now go!
        ser = TypeSerializer(class_name_map=dict(zip((dataset_to_name(ds) for ds in datasets), class_names)))
        for dataset, class_name in zip(datasets, class_names):
            ser.add_dataset(dataset)

        output = ser.emit(namespace=args.namespace, header_filename=os.path.basename(args.outfile_stub + ".h"))
        for suffix, output in output.items():
            files[suffix].write(output)
            print("Wrote file:", files[suffix].name)



