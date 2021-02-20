#ifndef GMRATREEIO_H
#define GMRATREEIO_H

#include "GMRANode.h"
#include "GMRAVisitor.h"
#include "GMRADecorator.h"
#include "GMRADataObject.h"
#include "NodeDistance.h"

#include <Eigen/Dense>

#include <list>
#include <vector>
#include <iostream>


template <typename TPrecision>
class GMRATreeIOWriter : public Visitor<TPrecision>{
  private:
    std::ostream &out;

  public:

    typedef typename GMRADataObject<TPrecision>::MatrixXp MatrixXp;
    typedef typename GMRADataObject<TPrecision>::VectorXp VectorXp;
    typedef typename GMRANode<TPrecision>::NodeVector NodeVector;


    GMRATreeIOWriter(std::ostream &outstream) : out(outstream){
    };

    //
    virtual void visit(GMRANode<TPrecision> *node){

      VectorXp center = node->getCenter();
      int size = center.size();
      out.write( (char*) &size, sizeof( int ) );
      out.write( (char*) center.data(), center.size()*sizeof( typename VectorXp::Scalar ) );

      int idim = node->getIntrinsicDimension();
      out.write( (char*) &idim, sizeof(int) );

      int scale = node->getScale();
      out.write( (char*) &scale, sizeof(int) );

      NodeVector &kids = node->getChildren();
      size = kids.size();
      out.write( (char*) &size, sizeof(int) );

      std::vector<int> &ind = node->getPoints();
      size = ind.size();
      out.write( (char*) &size, sizeof(int) );
      out.write( (char*) ind.data(), ind.size() * sizeof(int) );

      TPrecision radius = node->getRadius();
      out.write( (char*) &radius, sizeof(TPrecision) );

    };

};




//Tree data structure
template <typename TPrecision>
class GMRATreeIO{


public:
    typedef typename GMRADataObject<TPrecision>::MatrixXp MatrixXp;
    typedef typename GMRADataObject<TPrecision>::VectorXp VectorXp;




    static void writeTree( GMRATree<TPrecision> *tree, std::ostream &out){
      GMRATreeIOWriter< TPrecision > writer(out);
      tree->breadthFirstVisitor( &writer );
    };


    static GMRATree<TPrecision> *readTree( std::istream &in){

      IKMNode< TPrecision > *root = NULL;
      IKMNode< TPrecision > *current = NULL;
      std::list< IKMNode<TPrecision> *> nodes;
      std::list<int> kids;
      int nKidsCurrent = 0;

      while( in.peek() != EOF ){

        int dim;
        in.read( (char*) &dim, sizeof( int ) );
        VectorXp center( dim );
        in.read( (char*) center.data(), dim * sizeof( typename VectorXp::Scalar ) );

        int idim;
        in.read( (char*) &idim, sizeof(int) );


        int scale;
        in.read( (char*) &scale, sizeof(int) );

        int nKids;
        in.read( (char*) &nKids, sizeof(int) );

        int nPoints;
        in.read( (char*) &nPoints, sizeof(int) );
        std::vector<int> ind( nPoints );
        in.read( (char*) ind.data(), nPoints * sizeof(int) );

        TPrecision radius;
        in.read( (char*) &radius, sizeof(TPrecision) );

        IKMNode< TPrecision > *node = new IKMNode<TPrecision>(center, ind, radius, -1);
        node->setScale( scale );

        if( nKids > 0 ){
          nodes.push_front( node );
          kids.push_front( nKids );
        }

        if(current != NULL){
          current->addChild( node );
          node->setParent( current );
          nKidsCurrent--;
        }
        else{
          root = node;
        }

        if( nKidsCurrent == 0 && !nodes.empty() ){
          current = nodes.back();
          nodes.pop_back();
          nKidsCurrent = kids.back();
          kids.pop_back();
        }

      }

      IKMTree<TPrecision> *tree = new IKMTree<TPrecision>(NULL);
      tree->setRoot(root);
      return tree;

    };

};




#endif
