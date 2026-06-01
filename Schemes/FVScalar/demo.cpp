#include "FVmesh.hpp"

#include "mesh_builder.hpp"

#include <iostream>

/// The tools from FV are in their own namespace to avoid collisions (with terms such as Vertex/Edge...) 
// Most other tools are in the namespace Manicore, however it should only be useful when initializing the scheme
using namespace ManicoreFV;

/// The following meshes are already available
// The json files are not copied to build folder by make
// The relative path assumes that the program is run from a subfolder of the source tree
// The shared library is compiled with make, but the name assumes UNIX extension

/// Torus: flat torus. Uses 4 maps, the meshes are already periodic
// The only difference between the two shared libraries is the 3DEmbedding that may be used to visualize the solutions
// The first one embed into a plane (at z=0), while the second embed into a donught shaped torus.
// The change is only visual, the metric is still flat in any case.
[[maybe_unused]] constexpr const char * torusShared = "meshes/torus/libtorus_shared.so";
[[maybe_unused]] constexpr const char * torusShared3D = "meshes/torus/libtorus_3DEmbedding_shared.so";
// The following meshes are available. They use a cartesian grid
[[maybe_unused]] constexpr std::array torusMeshes{"../meshes/torus/torus_5.json",
                                                  "../meshes/torus/torus_10.json",
                                                  "../meshes/torus/torus_20.json",
                                                  "../meshes/torus/torus_25.json",
                                                  "../meshes/torus/torus_30.json",
                                                  "../meshes/torus/torus_35.json",
                                                  "../meshes/torus/torus_40.json"};
/// Sphere: surface of the unit sphere
// The meshes use 2 map: the north and south hemisphere, mapped to an unit disk with the stereographic projection
[[maybe_unused]] constexpr const char * sphereShared = "meshes/sphere/libsphere_shared.so";
[[maybe_unused]] constexpr std::array sphereMeshes{"../meshes/sphere/4_circle.json",
                                                   "../meshes/sphere/6_circle.json",
                                                   "../meshes/sphere/11_circle.json",
                                                   "../meshes/sphere/21_circle.json",
                                                   "../meshes/sphere/29_circle.json",
                                                   "../meshes/sphere/51_circle.json"};

int main() {
  // We first need to build the Manicore mesh object
  const char * meshfile = torusMeshes[0];
  const char * mapfile = torusShared;
  std::cout<<"Building mesh using json file \""<<meshfile<<"\", and map file: \""<<mapfile<<"\"\n";
  auto *geoMesh_p = Manicore::Mesh_builder<2>::build(meshfile,mapfile);
  // Then we can extract a FV mesh from the geoMesh
  FVMesh *mesh_p = FVMesh::build(geoMesh_p);
  // We no longer need the Manicore mesh (altought we may want to extract additionnal information for visualizing the solution).
  delete geoMesh_p;
  // FVMesh contains three lists: one of verticies, one of edges, and one of cells
  std::cout<<"The mesh contains "<<mesh_p->vertices.size()<<" vertices, "<<mesh_p->edges.size()<<" edges, "<<mesh_p->cells.size()<<" cells"<<std::endl;
  // The Vertex, Edge, and Cell objects all derive from baseCell
  Eigen::VectorXd uh = Eigen::VectorXd::Zero(mesh_p->cells.size());
  Eigen::VectorXd flux = Eigen::VectorXd::Ones(mesh_p->edges.size());
  // We can iterate over them using range-based for loop
  for (const Cell& T : mesh_p->cells) {
    const int iT = T.id(); // Global id, usable as index for uh
    const double cellArea = T.volume(); // Area of the cell (actual size on the manifold)
    double fl = 0.;
    for (size_t iET = 0; iET < T.edges().size(); ++iET) { // Loop over all edges of T
      const double wTE = T.boundaryOrientation()[iET]; // relative orientation of E with respect to the one of T
      const int iE = T.edges()[iET]->id(); // global index of the iET-th edge of T
      fl += wTE*flux(iE);
    }
    uh(iT) += fl/cellArea;
  }
  std::cout<<uh.transpose()<<std::endl;

  delete mesh_p;
  return 0;
}
