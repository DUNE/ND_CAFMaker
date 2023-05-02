import sys

import h5py
import numpy as np

numpy_to_hdf5 = {
	"i" : "H5T_STD_I",
	"l" : "H5T_STD_I",
	"B" : "H5T_STD_B",
	"f" : "H5T_IEEE_F",
	"d" : "H5T_IEEE_F",
}

class_template = """
struct {name}
{{
{members}
}};
"""

simple_member_template = "  {type} {name};"

enum_template = """
enum {name}_t
{{
{members}
}};

"""

enum_member_template = "  {name} = {val}, "

compound_type_template = """
H5::CompType BuildCompType<{klass}>()
{{
  H5::CompType ctype(sizeof({klass}));

{members}

  return ctype;
}}
"""

compound_type_member_template = '  ctype.insertMember("{name}", HOFFSET({klass}, {name}), {h5type});'

class TypeSerializer:
	def __init__(self, ds):
		self.dataset = ds

		self.discovered_enums = {}
		self.ctype_members = []

	def EnumString(self):
		ret = ""
		for name, enum in self.discovered_enums.items():
			members = []
			for itemname, val in enum.items():
				members.append(enum_member_template.format(name=itemname, val=val))
			ret += enum_template.format(name=name, members="\n".join(members))
		return ret

	def TypeString(self, typ, fieldname=None, which="cpp"):
		assert which in ("cpp", "h5"), "Unrecognized kind of type string: " + which

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
			h5_name = "H5::PredType::" + self.TypeString(typ, fieldname=None)
		elif h5py.check_enum_dtype(typ):
			typenames = fieldname.split("_")
			typename = "".join(t.capitalize() for t in typenames)
			self.discovered_enums[typename] = { "k{name}".format(name="".join(k.split())): v for k, v in h5py.check_enum_dtype(typ).items() }
			h5_name = "H5::EnumType(H5::PredType::{typ})".format(typ=self.TypeString(np.dtype(str(typ)), fieldname))
			# print("discovered enum:", typename, self.discovered_enums[typename])
			cpp_name += typename
		elif h5py.check_vlen_dtype(typ):
			# for variable-length items, we have to maintain a special "handle"
			base_type = self.TypeString(h5py.check_vlen_dtype(typ), fieldname)
			cpp_name += "hvl_t {name}_handle;  std::vector<{base_type}>".format(name=fieldname, base_type=base_type)
			h5_name = "H5::VarLenType(H5::PredType::{typ})".format(typ=base_type)
		elif typ.ndim > 0:
			assert typ.ndim == 1, "Don't know how to handle multi-dimensional types"
			base_type = self.TypeString(typ.subdtype[0], fieldname)
			cpp_name += "std::array<{type}, {len}>".format(type=base_type, len=typ.shape[0])
			h5_name = "H5::ArrayType(H5::PredType::{typ}, 1, &std::array<hsize_t, 1>{{{count}}}[0])".format(typ=base_type, count=typ.shape[0])
		elif typ.char in numpy_to_hdf5:
			cpp_name += numpy_to_hdf5[typ.char]
			cpp_name += str(typ.itemsize * 8)
			cpp_name += "BE" if typ.byteorder == ">" else "LE" if typ.byteorder == "<" else sysorder[sys.byteorder]
			h5_name = "H5::PredType::" + cpp_name
		else:
			raise TypeError("Don't know how to handle type: " + str(typ))

		h5_name = h5_name.replace("H5T_", "")

		return cpp_name if which == "cpp" else h5_name

	def Serialize(self, typename):
		ret = ""

		members = []
		h5_members = []
		for fieldname in self.dataset.dtype.names:
			typ = self.dataset.dtype[fieldname]
			members.append(simple_member_template.format(name=fieldname, type=self.TypeString(typ, fieldname=fieldname)))
			h5_members.append(compound_type_member_template.format(name=fieldname, klass=typename, h5type=self.TypeString(typ, fieldname=fieldname, which="h5")))

		ret += self.EnumString()

		ret += class_template.format(name=typename, members="\n".join(members))
		ret += compound_type_template.format(klass=typename, members="\n".join(h5_members))

		return ret
