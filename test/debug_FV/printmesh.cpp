#include "FVmesh.hpp"

#include "mesh_builder.hpp"
#include "quadraturerule.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace ManicoreFV;

// Include a constraint to resolve ambiguous overload
template<typename Derived>
std::ostream & operator<<(std::ostream &out, const Eigen::DenseBase<Derived>& v) requires (static_cast<bool>(Derived::IsVectorAtCompileTime)) {
  out<<v(0);
  for (int i = 1; i < v.size(); ++i) {
    out<<", "<<v(i);
  }
  return out;
}
template<typename T1>
std::ostream & operator<<(std::ostream &out, const std::vector<T1*>& v) {
  out<<v[0]->id();
  for (size_t i = 1; i < v.size(); ++i) {
    out<<", "<<v[i]->id();
  }
  return out;
}
template<typename T1>
std::ostream & operator<<(std::ostream &out, const std::vector<T1>& v) {
  out<<v[0];
  for (size_t i = 1; i < v.size(); ++i) {
    out<<", "<<v[i];
  }
  return out;
}
void printMesh(const char * meshfile, const char * mapfile, const char *prefix) {
  std::cout<<meshfile<<"\n";
  auto *geoMesh_p = Manicore::Mesh_builder<2>::build(meshfile,mapfile);
  FVMesh *mesh_p = FVMesh::build(geoMesh_p);
  std::cout<<"Mesh with: "<<mesh_p->vertices.size()<<" vertices, "<<mesh_p->edges.size()<<" edges, "<<mesh_p->cells.size()<<" cells"<<std::endl;
  // Print midpoints
  {
    // MapID, LocX, LocY, Id
    std::ofstream fh("vert.csv");
    for (const Vertex& V : mesh_p->vertices) {
      const Manicore::dCell_map VGeo = geoMesh_p->get_cell_map<0>(V.id());
      for (size_t rId = 0; rId < geoMesh_p->get_map_ids(0,V.id()).size(); ++rId) {
        fh<<geoMesh_p->get_map_ids(0,V.id())[rId]<<", "<<VGeo.coord(rId)<<", "<<V.id()<<"\n";
      }
    }
  }
  auto computeCenter = [geoMesh_p](const auto& T,size_t mId)->Eigen::Vector2d {
    Eigen::Vector2d rv(0.,0.);
    for (Vertex *V : T.vertices()) {
      const Manicore::dCell_map VGeo = geoMesh_p->get_cell_map<0>(V->id());
      const std::vector<size_t> &vmids = geoMesh_p->get_map_ids(0,V->id());
      for (size_t i = 0; i < vmids.size(); ++i) {
        if (vmids[i] == mId) {
          rv += VGeo.coord(i);
          break;
        }
        if (i == vmids.size()-1) {
          std::cerr<<"Warning: Failed to find corresponding map Id in vertex "<<V->id()<<" of cell "<<T.id()<<std::endl;
          std::cerr<<"Searching for map Id "<<mId<<", but the vertex maps are "<<vmids<<std::endl;
        }
      }
    }
    return rv/T.vertices().size();
  };
  auto getMaps = [geoMesh_p]<typename T>(const T& e)->std::vector<size_t> {
    if constexpr (std::is_same<T,Vertex>::value) {
      return geoMesh_p->get_map_ids(0,e.id());
    } else if constexpr (std::is_same<T,Edge>::value) {
      return geoMesh_p->get_map_ids(1,e.id());
    } else {
      return geoMesh_p->get_map_ids(2,e.id());
    }
  };
  {
    // MapID, LocX, LocY, Id
    std::ofstream fh("edge.csv");
    for (const Edge& E : mesh_p->edges) {
      for (size_t mId : getMaps(E)) {
        fh<<mId<<", "<<computeCenter(E,mId)<<", "<<E.id()<<"\n";
      }
    }
  }
  {
    // MapID, LocX, LocY, Id
    std::ofstream fh("cell.csv");
    for (const Cell& T : mesh_p->cells) {
      for (size_t mId : getMaps(T)) {
        fh<<mId<<", "<<computeCenter(T,mId)<<", "<<T.id()<<"\n";
      }
    }
  }
  {
    // EdgeId, V1, V2, O1, O2
    std::ofstream fh("edgeVertices.csv");
    for (const Edge& E : mesh_p->edges) {
      fh<<E.id()<<", "<<E.vertices()<<", "<<E.boundaryOrientation()<<"\n";
    }
  }
  {
    // CellId, isFlat, E1, E2, ... , O1, O2, ...
    std::ofstream fh("cellEdges.csv");
    for (const Cell& T : mesh_p->cells) {
      const auto &TGeo = geoMesh_p->get_cell_map<2>(T.id());
      fh<<T.id()<<", "<<(TGeo.is_flat()? 1 : 0)<<", "<<T.edges()<<", "<<T.boundaryOrientation()<<"\n";
    }
  }
  {
    // MapId, CellId, LocX1, LocY1, LocX2, LocY2, ...
    std::ofstream fh("quadlocs.csv");
    // MapId, CellId, Basis1Loc1, Basis2Loc1, ...
    std::ofstream fhPoly("polyeval.csv"); 
    constexpr int degQuad = 5;
    constexpr int degPoly = 2;
    Manicore::Monomial_powers<2>::init(degPoly);
    {
      std::ofstream fhBasis("polyname.csv");
      fhBasis << "1";
      std::vector<std::array<size_t,2>> basis = Manicore::Monomial_powers<2>::complete(degPoly);
      for (size_t pBasis = 1; pBasis< basis.size();++pBasis){
        fhBasis << ", ";
        if (basis[pBasis][0] > 0) {
          fhBasis<<"X^"<<basis[pBasis][0];
        }
        if (basis[pBasis][1] > 0) {
          fhBasis<<"Y^"<<basis[pBasis][1];
        }
      }
      fhBasis<<"\n";
    }
    for (const Cell& T : mesh_p->cells) {
      const auto &TGeo = geoMesh_p->get_cell_map<2>(T.id());
      Manicore::QuadratureRule<2> quad = Manicore::generate_quadrature_rule(TGeo, degQuad);
      for (size_t itMId = 0; itMId < geoMesh_p->get_map_ids(2,T.id()).size(); ++itMId) {
        const size_t mId = geoMesh_p->get_map_ids(2,T.id())[itMId];
        fh<<mId<<", "<<T.id();
        fhPoly<<mId<<", "<<T.id();
        for (const auto& node : quad) {
          Eigen::Vector2d x = TGeo.evaluate_I(itMId, node.vector);
          fh<<", "<<x;
          for (size_t pBasis = 0; pBasis < Manicore::Dimension::PolyDim(degPoly,2); ++pBasis) {
            fhPoly<<", "<<TGeo.evaluate_poly_pullback(itMId,x,pBasis,degPoly);
          }
        }
      }
      fh<<"\n";
      fhPoly<<"\n";
    }
  }
  { // Convert to image
    std::system((std::string("python3 printmesh.py ")+prefix).c_str());
    std::remove("vert.csv");
    std::remove("edge.csv");
    std::remove("cell.csv");
    std::remove("edgeVertices.csv");
    std::remove("cellEdges.csv");
    std::remove("quadlocs.csv");
    std::remove("polyeval.csv");
    std::remove("polyname.csv");
  }
  { // Print area
    double totalArea = 0.;
    for (const Cell& T : mesh_p->cells) {
      std::cout<<"Cell: "<<T.id()<<", volume: "<<T.volume()<<"\n";
      totalArea += T.volume();
    }
    std::cout<<"Total area: "<<totalArea<<std::endl;
  }
  
  delete geoMesh_p;
  delete mesh_p;
}
 int main() {
  printMesh("../meshes/test/23_pts.json","meshes/test/libdisk_maps.so","test23");
  printMesh("../meshes/torus/torus_5.json","meshes/torus/libtorus_3DEmbedding_shared.so","T5");
  printMesh("../meshes/torus/torus_40.json","meshes/torus/libtorus_3DEmbedding_shared.so","T40");
  printMesh("../meshes/sphere/4_circle.json","meshes/sphere/libsphere_shared.so","S4");
  printMesh("../meshes/sphere/29_circle.json","meshes/sphere/libsphere_shared.so","S29");
  return 0;
}
