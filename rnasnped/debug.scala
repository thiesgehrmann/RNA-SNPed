import scala.tools.nsc.interpreter.{ILoop, SimpleReader}
import scala.tools.nsc.Settings



object debug {

  def breakIf( test: Boolean) = {
    val repl = new ILoop
    repl.settings = new Settings
    repl.in = SimpleReader()
    
    // set the "-Yrepl-sync" option
    repl.settings.Yreplsync.value = true
    
    // start the interpreter and then close it after you :quit
    repl.createInterpreter()
    repl.loop()
    repl.closeInterpreter()
  }

}
