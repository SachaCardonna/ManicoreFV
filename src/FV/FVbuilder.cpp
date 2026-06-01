#include "FVmesh.hpp"

#include "mesh.hpp"
#include "integral.hpp"

using namespace ManicoreFV;

FVMesh* FVMesh::build(const Manicore::Mesh<2> *mesh) {
  FVMesh *fvm = new FVMesh();
  const size_t nbV = mesh->n_cells(0);
  const size_t nbE = mesh->n_cells(1);
  const size_t nbT = mesh->n_cells(2);
  // Retrieve connectivity from Mesh
  std::vector<Vertex*> vVert(nbV);
  std::vector<Edge*> vEdge(nbE);
  for (size_t iV = 0; iV < nbV; ++iV) {
    fvm->vertices.push_back(iV);
    fvm->vertices.back()._vol = 1.;
    vVert[iV] = &(fvm->vertices.back());
  }
  for (size_t iE = 0; iE < nbE; ++iE) {
    fvm->edges.push_back(iE);
    vEdge[iE] = &(fvm->edges.back());
    for (size_t iVE = 0; iVE < mesh->get_boundary(0,1,iE).size(); iVE++) {
      const size_t iVEg = mesh->get_boundary(0,1,iE)[iVE];
      fvm->edges.back()._V.push_back(vVert[iVEg]);
      vVert[iVEg]->_V.push_back(&(fvm->edges.back()));
      fvm->edges.back()._Borientation.push_back(mesh->get_boundary_orientation(1,iE,iVE));
    }
  }
  for (size_t iT = 0; iT < nbT; ++iT) {
    fvm->cells.push_back(iT);
    for (size_t iVTg : mesh->get_boundary(0,2,iT)) {
      fvm->cells.back()._V.push_back(vVert[iVTg]);
      vVert[iVTg]->_E.push_back(&(fvm->cells.back()));
    }
    for (size_t iET = 0; iET < mesh->get_boundary(1,2,iT).size(); iET++) {
      const size_t iETg = mesh->get_boundary(1,2,iT)[iET];
      fvm->cells.back()._E.push_back(vEdge[iETg]);
      vEdge[iETg]->_E.push_back(&(fvm->cells.back()));
      fvm->cells.back()._Borientation.push_back(mesh->get_boundary_orientation(2,iT,iET));
    }
  }
  // Compute area
  constexpr int DOE = 5;
  for (Edge& E : fvm->edges) {
    Manicore::Integral<2,1> integral(mesh);
    const auto quad = integral.generate_quad(E.id(),DOE);
    const Eigen::VectorXd vol = integral.evaluate_volume_form(E.id(),quad);
    double area = 0.;
    for (size_t iqn = 0; iqn < quad.size(); ++iqn) {
      area += vol(iqn)*quad[iqn].w;
    }
    E._vol = area;
  }
  for (Cell& T : fvm->cells) {
    Manicore::Integral<2,2> integral(mesh);
    const auto quad = integral.generate_quad(T.id(),DOE);
    const Eigen::VectorXd vol = integral.evaluate_volume_form(T.id(),quad);
    double area = 0.;
    for (size_t iqn = 0; iqn < quad.size(); ++iqn) {
      area += vol(iqn)*quad[iqn].w;
    }
    T._vol = area;
  }
  return fvm;
};
