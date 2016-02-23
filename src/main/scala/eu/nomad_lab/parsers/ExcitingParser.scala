package eu.nomad_lab.parsers
import eu.nomad_lab.DefaultPythonInterpreter
import org.{json4s => jn}

object ExcitingParser extends SimpleExternalParserGenerator(
  name = "ExcitingParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("ExcitingParser")) ::
      ("version" -> jn.JString("1.0")) :: Nil),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """\s*=================================================+\s*
\s*\|\s*EXCITING\s(?<version>\S*) started\s*=
(?:\s*\|\sversion hash id:\s*(?<hashId>\S*)\s*=)?""".r,
  cmd = Seq(DefaultPythonInterpreter.python2Exe(), "${envDir}/parsers/exciting/parser/parser-exciting/parser_exciting.py",
    "--uri", "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-exciting/parser_exciting.py",
    "parser-exciting/setup_paths.py",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/exciting.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-exciting" -> "parsers/exciting/parser/parser-exciting/",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info",
    "python" -> "python-common/common/python/nomadcore") ++ DefaultPythonInterpreter.commonDirMapping()
)
