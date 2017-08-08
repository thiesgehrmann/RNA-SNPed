import scala.io.Source
import java.io._

///////////////////////////////////////////////////////////////////////////////
package rnasnped {
object sampleTree {

  case class treeNode(id: Int = -1, name: String = "", children: Array[Int] = Array.empty[Int], isLeaf: Boolean = false, height: Int = -1) {

  }

  /////////////////////////////////////////////////////////////////////////////

  def readTree(treeFile: String) = {

    case class treeNodeParserLine(id: Int, name: String, children: Array[Int], isLeaf: Boolean)

    def parseTreeFileLine(line: String) = {
      val s        = line.split('\t');
      val id       = s(0).toInt - 1 //to offset the 1 index in the file
      val name     = s(1)
      val children = if(s(2) == "-") Array.empty[Int] else s(2).split(',').map(x => x.toInt -1);
      treeNodeParserLine(id, name, children, s(2) == "-")
    }
    
    def getHeight(tree: Array[treeNodeParserLine], node: Int): Array[treeNode] ={
      val treeNodeParserLine(id, name, children, isLeaf) = tree(node);
      if(isLeaf){
        Array(treeNode(id, name, children, true, 1))
      } else {
        val childNodes = children.map( x => getHeight(tree, x)).flatten
        childNodes :+ treeNode(id, name, children, false, childNodes.last.height + 1)
      }
    }
    
    val parsedLines = io.Source.fromFile(treeFile).getLines.filter(x => x.length > 0 && x(0) != '#').map( x => parseTreeFileLine(x)).toArray

    getHeight(parsedLines, parsedLines.length - 1).sortBy(_.id)
  }

  /////////////////////////////////////////////////////////////////////////////

  def parentPaths(tree: Array[treeNode]) = {
    val parents = collection.mutable.ArrayBuffer.fill(tree.length){-1}
    tree.zipWithIndex.foreach{ case (node, node_i) =>
      node.children.foreach{ child_i => parents(child_i) = node_i}
    }

    def parentPathHelper(nodeID: Int, parentsList: List[Int]): List[Int] = {
      parents(nodeID) match {
        case -1 => parentsList :+ nodeID
        case _  => parentPathHelper(parents(nodeID), parentsList :+ parents(nodeID) :+ nodeID)
      }
    }

    tree.indices.map( node_i => parentPathHelper(node_i, List.empty[Int]))
  }

  /////////////////////////////////////////////////////////////////////////////

  def getNodeLeaves(tree: Array[treeNode]) = {

    def getNodeLeavesHelper(path: Array[Int], visited: Set[Int], leaves: Array[Array[Int]]): Array[Array[Int]] = {
      if (path.length > 0) {
        val nodeID = path(0)
        val unvisitedChildren = tree(nodeID).children.filter(c => !visited.contains(c))

        val newLeaves = if (tree(nodeID).isLeaf) {
            leaves.zipWithIndex.map{ case (nl, nid) => if (path contains nid) nl :+ nodeID else nl }
          } else {
            leaves
          }

        unvisitedChildren.length match {
          case 0 => {
            val pathTail = path.tail
            getNodeLeavesHelper(pathTail, visited + nodeID, newLeaves)
          }
          case _ => {
            getNodeLeavesHelper(unvisitedChildren.head +: path, visited + nodeID, leaves)
          }
        }
 
      } else {
        leaves
      }
      
    }

    getNodeLeavesHelper(Array(tree.length-1), Set.empty[Int], tree.map(n => Array.empty[Int]))

  }

  /////////////////////////////////////////////////////////////////////////////

}
}
