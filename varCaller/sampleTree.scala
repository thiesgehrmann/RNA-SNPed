import scala.io.Source
import java.io._

///////////////////////////////////////////////////////////////////////////////
package rnasnped {
object sampleTree {

  case class treeNode(id: Int = -1, name: String = "", children: Array[Int] = Array.empty[Int], isLeaf: Boolean = false, height: Int = -1)

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

}
}
