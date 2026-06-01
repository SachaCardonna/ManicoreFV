#ifndef FVMESH_HPP_INCLUDED
#define FVMESH_HPP_INCLUDED

#include <vector>
#include <list>
#include <Eigen/Dense>

namespace Manicore {
  template<size_t> class Mesh;
}

/** @defgroup FiniteVolume
  @brief Classes providing support for finite volumes schemes
  */

namespace ManicoreFV {
  
  /** \file FVmesh.hpp
    Generate connectivity graph and specialized integration rule on edges
    */

  class Vertex;
  class Edge;
  class Cell;
  struct FVMesh;

  /// \addtogroup FiniteVolume
  ///@{

  /// Base class for all cell types
  /** The connectivity information is stored using pointer to the boundary elements
    */
  template<typename T1, typename T2>
  class baseCell {
    public:
      /// Return the global id of the given cell
      int id() const {return _Id;}
      const std::vector<Vertex*> &vertices() const requires(std::same_as<T1,Vertex>) {return _V;}
      const std::vector<Edge*> &edges() const requires(std::same_as<T1,Edge> || std::same_as<T2,Edge>) {if constexpr (std::same_as<T1,Edge>) return _V; else return _E;}
      const std::vector<Cell*> &cells() const requires(std::same_as<T2,Cell>) {return _E;}
      /// Relative orientation of boundary
      const std::vector<float>& boundaryOrientation() const requires(std::same_as<T1,Vertex>) {return _Borientation;}
      /// Volume of the cell
      double volume() const {return _vol;}
      // TODO implement an integrator
    private:
      baseCell(int id) : _Id(id) {}
      int _Id;
      std::vector<T1*> _V;
      std::vector<T2*> _E;
      std::vector<float> _Borientation;
      double _vol; // Volume of the cell, TODO replace with integrator when doing high-order FV
      friend FVMesh;
  };

  /// Specialization for vertices
  class Vertex : public baseCell<Edge,Cell> {
    using baseCell<Edge,Cell>::baseCell;
  };

  /// Specialization for edges
  class Edge : public baseCell<Vertex,Cell> {
    using baseCell<Vertex,Cell>::baseCell;
  };

  /// Specialization for cells
  class Cell : public baseCell<Vertex,Edge> {
    using baseCell<Vertex,Edge>::baseCell;
  };

  /// Structure storing the mesh
  struct FVMesh {
    std::list<Vertex> vertices;
    std::list<Edge> edges;
    std::list<Cell> cells;
    static FVMesh* build(const Manicore::Mesh<2>*);
  };

  ///@}

}

#endif

