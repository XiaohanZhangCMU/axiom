/*
 * Defines python interfaces to C++ classes and methods
 */

#include <cstdio>
#include <iostream>
#include <pybind11/pybind11.h>

#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include "fem.h"

PYBIND11_MAKE_OPAQUE(std::vector<int>);

namespace bp = pybind11;

std::string FEMFrame::get_contact_file() { std::string str(contact_file); return str; }
void FEMFrame::set_contact_file(std::string s) { strcpy(contact_file, s.c_str()); }

std::vector<int> FEMFrame::get_bNds() {return std::vector<int>(bNds,bNds+n_bdy_nodes);}
void FEMFrame::set_bNds(std::vector<int> v) {
    if (bNds==NULL)  { 
        bNds = (int*) realloc(bNds, sizeof(int)*v.size());
        n_bdy_nodes = v.size();
    } 
    memcpy(bNds, &v[0], sizeof(int)*n_bdy_nodes); 
}

std::vector<int> FEMFrame::get_bTags() {return std::vector<int>(bTags,bTags+n_bdy_nodes);}
void FEMFrame::set_bTags(std::vector<int> v) {
    if (bTags==NULL)  { 
        bTags = (int*) realloc(bTags, sizeof(int)*v.size());
        n_bdy_nodes = v.size();
    } 
    memcpy(bTags, &v[0], sizeof(int)*n_bdy_nodes); 
}

class Py_SCParser : public SCParser { 
    using SCParser::SCParser;
    int assignvar(int offset) override { PYBIND11_OVERLOAD(int, SCParser, assignvar, offset); } 
};

class Py_Organizer : public Organizer { 
    using Organizer::Organizer;
};

class Py_MDFrame : public MDFrame {
public:
    using MDFrame::MDFrame;
    void Alloc() override { PYBIND11_OVERLOAD(void, MDFrame, Alloc); }
    void call_potential() override { PYBIND11_OVERLOAD(void, MDFrame,call_potential); } 
    void potential() override { PYBIND11_OVERLOAD(void, MDFrame, potential); } 
    void potential_energyonly() override { PYBIND11_OVERLOAD(void, MDFrame, potential_energyonly); } 
    double potential_energyonly(int iatom) override { PYBIND11_OVERLOAD(double, MDFrame, potential_energyonly, iatom); } 
    void plot() override { PYBIND11_OVERLOAD(void, MDFrame, plot); }
}; 

class Py_FEMFrame : public FEMFrame {
public:
    using FEMFrame::FEMFrame;
    void Alloc() override { PYBIND11_OVERLOAD(void, FEMFrame, Alloc); }
    void potential() override { PYBIND11_OVERLOAD(void, FEMFrame, potential); } 
    void plot() override { PYBIND11_OVERLOAD(void, MDFrame, plot); }
}; 

namespace axiom {

PYBIND11_MODULE(mdfem, m) {
    bp::bind_vector<std::vector<int>>(m, "VectorInt", bp::buffer_protocol());
    bp::bind_vector<std::vector<double>>(m, "VectorDouble", bp::buffer_protocol());

    // optional module doc string
    m.doc() = "Tensor class exposure";

    bp::class_<SCParser, Py_SCParser>(m, "SCParser")
        .def(bp::init<> ())
        ;

    bp::class_<Organizer, SCParser, Py_Organizer>(m, "Organizer")
        .def(bp::init<> ())
        .def("setnolog", &Organizer::setnolog)
        .def("sleep", &Organizer::getsleep)
        .def_property("dirname", &Organizer::get_dirname, &Organizer::set_dirname, "property: dirname")
        ;

    bp::class_<MDFrame, Organizer, Py_MDFrame>(m, "MDFrame")
        .def(bp::init<> ())
        .def_property("incnfile", &MDFrame::get_incnfile, &MDFrame::set_incnfile, "property: incnfile")
        .def_property("finalcnfile", &MDFrame::get_finalcnfile, &MDFrame::set_finalcnfile, "property: finalcnfile")
        .def_property("crystalstructure", &MDFrame::get_crystalstructure, &MDFrame::set_crystalstructure, "property: crystalstructure")
        .def_property("latticeconst", &MDFrame::get_latticeconst, &MDFrame::set_latticeconst, "property: latticeconst")
        .def_property("latticesize", &MDFrame::get_latticesize, &MDFrame::set_latticesize, "property: latticesize")

        .def_property("atomradius", &MDFrame::get_atomradius, &MDFrame::set_atomradius, "property: atomradius")
        .def_property("atomcolor", &MDFrame::get_atomcolor, &MDFrame::set_atomcolor, "property: atomcolor")
        .def_property("bondcolor", &MDFrame::get_bondcolor, &MDFrame::set_bondcolor, "property: bondcolor")
        .def_property("backgroundcolor", &MDFrame::get_backgroundcolor, &MDFrame::set_backgroundcolor, "property: backgroundcolor")
        .def_property("fixatomcolor", &MDFrame::get_fixatomcolor, &MDFrame::set_fixatomcolor, "property: fixatomcolor")
        .def_property("highlightcolor", &MDFrame::get_highlightcolor, &MDFrame::set_highlightcolor, "property: hightlightcolor")
        .def_property("rotateangles", &MDFrame::get_rotateangles, &MDFrame::set_rotateangles, "property: roateangles")
        .def_property("plot_color_windows", &MDFrame::get_plot_color_windows, &MDFrame::set_plot_color_windows, "property: plot_color_windows")
        .def_property("plot_limits", &MDFrame::get_plot_limits, &MDFrame::set_plot_limits, "property: plot_limits")

        .def("call_potential", &MDFrame::call_potential)
        .def("eval", &MDFrame::eval)
        .def("makecrystal", &MDFrame::makecrystal)
        .def("writecn", &MDFrame::writefinalcnfile)
        .def("openwin", &MDFrame::openwin)
        .def("relax",   &MDFrame::relax)
        .def("plot",    &MDFrame::plot)
        .def("saverot", &MDFrame::saverot)
        .def("rotate",  &MDFrame::rotate)
        .def("alloccolors", &MDFrame::alloccolors)
        .def("refreshnnlist", &MDFrame::NbrList_refresh)

        .def_readwrite("plotfreq", &MDFrame::plotfreq)
        .def_readwrite("plot_atom_info", &MDFrame::plot_atom_info)
        .def_readwrite("bondradius", &MDFrame::bondradius)
        .def_readwrite("bondlength", &MDFrame::bondlength)
        .def_readwrite("conj_ftol",  &MDFrame::conj_ftol)
        .def_readwrite("conj_itmax", &MDFrame::conj_itmax)
        .def_readwrite("conj_fevalmax", &MDFrame::conj_fevalmax)
        .def_readwrite("conj_fixbox",   &MDFrame::conj_fixbox)
        ;

    bp::class_<FEMFrame, MDFrame, Py_FEMFrame>(m,"FEMFrame", bp::buffer_protocol())
        .def(bp::init<> ())
        .def_readwrite("n_bdy_nodes", &FEMFrame::n_bdy_nodes)
        .def_property("contact_file", &FEMFrame::get_contact_file, &FEMFrame::set_contact_file, "property: contact_file")
        .def("get_contact_file", &FEMFrame::get_contact_file)
        .def_property("bNds", &FEMFrame::get_bNds, &FEMFrame::set_bNds, "property: bNds")
        .def_property("bTags", &FEMFrame::get_bTags, &FEMFrame::set_bTags, "property: bTags")
        ;

} /* PYBIND11_MODULE */

} /* namespace axiom */

