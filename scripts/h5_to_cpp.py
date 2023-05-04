import argparse
import contextlib
import itertools
import os.path
import sys

import h5py
import numpy as np

# -------------------------------------------------------

numpy_to_hdf5 = {
    "i": "H5T_STD_I",
    "l": "H5T_STD_I",
    "B": "H5T_STD_B",
    "f": "H5T_IEEE_F",
    "d": "H5T_IEEE_F",
}

# -------------------------------------------------------

hdr_template = \
"""
//  WARNING:
//    This file was autogenerated by {filename}.
//    Do not edit by hand!

#include <array>
#include <vector>

#include "H5Cpp.h"

namespace {namespace}
{{

  // This generic template will be overloaded
  // for every specific type that we create below
  template <typename T>
  H5::CompType BuildCompType();

{members}

}}
"""

# -----------

impl_template = \
"""
//  WARNING:
//    This file was autogenerated by {filename}.
//    Do not edit by hand!

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

enum_template = \
"""
enum {name}
{{
{members}
}};
"""

# -----------

enum_member_template = "{name} = {val},"

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

compound_type_member_template = 'ctype.insertMember("{name}", HOFFSET({klass}, {name}), {h5type});'


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
    def __init__(self):
        self.discovered_enums = {}
        self.cpp_types = {}
        self.comptype_builders_decl = {}
        self.comptype_builders_impl = {}

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

        if h5py.check_dtype(ref=typ):
            cpp_name += "hdset_reg_ref_t"
            h5_name = "H5::PredType::STD_REF_DSETREG"
        # horrible hack, but bools are always stored as ints, so we have no other way of knowing
        elif fieldname and fieldname.startswith("is_"):
            cpp_name += "bool"
            h5_name = "H5::PredType::" + self.type_string(typ, fieldname=None)
        elif h5py.check_enum_dtype(typ):
            typenames = fieldname.split("_")
            typename = "".join(t.capitalize() for t in typenames) + "_t"

            # we really shouldn't have the side effect of storing the enums
            # performed inside a function that's doing something else,
            # but it's the simplest way to do it
            self.discovered_enums[typename] = Serializable(template=enum_template, template_args={"name": typename},
                                                           member_list=[Serializable(template=enum_member_template,
                                                                                     template_args=dict(
                                                                                         name="k{name}".format(
                                                                                             name="".join(k.split())),
                                                                                         val=v))
                                                                        for k, v in h5py.check_enum_dtype(typ).items()])

            h5_name = "H5::EnumType(H5::PredType::{typ})".format(typ=self.type_string(np.dtype(str(typ)), fieldname))
            # print("discovered enum:", typename, self.discovered_enums[typename])
            cpp_name += typename
        elif h5py.check_vlen_dtype(typ):
            base_type = self.type_string(h5py.check_vlen_dtype(typ), fieldname)
            cpp_name += "std::vector<{base_type}>".format(base_type=base_type)
            h5_name = "H5::VarLenType(H5::PredType::{typ})".format(typ=base_type)
        elif typ.ndim > 0:
            assert typ.ndim == 1, "Don't know how to handle multi-dimensional types"
            base_type = self.type_string(typ.subdtype[0], fieldname)
            cpp_name += "std::array<{type}, {len}>".format(type=base_type, len=typ.shape[0])
            h5_name = "H5::ArrayType(H5::PredType::{typ}, 1, &std::array<hsize_t, 1>{{{count}}}[0])".format(
                typ=base_type, count=typ.shape[0])
        elif typ.char in numpy_to_hdf5:
            cpp_name += numpy_to_hdf5[typ.char]
            cpp_name += str(typ.itemsize * 8)
            cpp_name += "BE" if typ.byteorder == ">" else "LE" if typ.byteorder == "<" else sysorder[sys.byteorder]
            h5_name = "H5::PredType::" + cpp_name
        else:
            raise TypeError("Don't know how to handle type: " + str(typ))

        h5_name = h5_name.replace("H5T_", "")

        return cpp_name if which == "cpp" else h5_name

    def add_dataset(self, dataset, class_name):
        self._dirty = True

        cpp_members = []
        cpp_private_members = []
        h5_members = []
        print("Examining dataset:", dataset.name)
        for fieldname in dataset.dtype.names:
            typ = dataset.dtype[fieldname]
            cpp_members.append(Serializable(template=simple_member_template,
                                            template_args=dict(name=fieldname,
                                                               typ=self.type_string(typ,fieldname=fieldname))))

            # for variable-length items, we have to maintain a special "handle"
            # in addition to the vector generated for cpp_members
            if h5py.check_vlen_dtype(typ):
                cpp_private_members.append(Serializable(template=simple_member_template,
                                                        template_args=dict(name=fieldname + "_handle", typ="hvl_t"),
                                                        base_indent="  "))

            h5_members.append(Serializable(template=compound_type_member_template,
                                           template_args=dict(name=fieldname,
                                                              klass=class_name,
                                                              h5type=self.type_string(typ, fieldname=fieldname, which="h5"))))
        if len(cpp_private_members) > 0:
            cpp_members.append(Serializable(template="\nprivate:", template_args={}))
            cpp_members += cpp_private_members

        # todo: need to add `template <> const hdset_reg_ref_t& GetRef<Particle>() const { return particles; }` (etc.)
        #       but ONLY if they're structured datatypes...

        self.cpp_types[class_name] = Serializable(template=class_template, template_args=dict(name=class_name),
                                                  member_list=cpp_members)
        self.comptype_builders_decl[class_name] = Serializable(template=compound_type_decl_template,
                                                               template_args=dict(klass=class_name))
        self.comptype_builders_impl[class_name] = Serializable(template=compound_type_impl_template,
                                                               template_args=dict(klass=class_name),
                                                               member_list=h5_members)

    def emit(self, namespace, header_filename):
        if len(self._serializables) == 0 or self._dirty:
            self._serializables[".h"] = Serializable(template=hdr_template,
                                                     template_args=dict(filename=os.path.basename(__file__), namespace=namespace),
                                                     member_list=tuple(itertools.chain(self.discovered_enums.values(),
                                                                                       self.cpp_types.values(),
                                                                                       self.comptype_builders_decl.values())),
                                                     member_indent="  ")
            self._serializables[".cxx"] = Serializable(template=impl_template,
                                                       template_args=dict(filename=os.path.basename(__file__), namespace=namespace, hdr=header_filename),
                                                       member_list=self.comptype_builders_impl.values(),
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

    # todo: need to ensure `events` dataset always comes last.
    #       otherwise we'd need to have forward declarations for the other types
    #       so that the GetRef<>() methods that are emitted don't refer to unknown types.
    #       (could try to do more generally and just put types that have region refs to datasets of structured types last,
    #        but that sounds too hard right now)

    with contextlib.ExitStack() as stack:
        # this opens all of the files with a context manager
        # (so they all get closed properly, even if there's an exception)
        files = {suffix: stack.enter_context(open(args.outfile_stub + suffix, "w")) for suffix in (".h", ".cxx")}

        # now go!
        ser = TypeSerializer()
        for dataset, class_name in zip(datasets, class_names):
            ser.add_dataset(dataset, class_name=class_name)

        output = ser.emit(namespace=args.namespace, header_filename=os.path.basename(args.outfile_stub + ".h"))
        for suffix, output in output.items():
            files[suffix].write(output)
            print("Wrote file:", files[suffix].name)



