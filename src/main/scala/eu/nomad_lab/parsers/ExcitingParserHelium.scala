package eu.nomad_lab.parsers

import eu.{ nomad_lab => lab }
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import scala.collection.breakOut

object ExcitingParserHelium extends SimpleExternalParserGenerator(
  name = "ExcitingParserHelium",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("ExcitingParserHelium")) ::
      ("parserId" -> jn.JString("ExcitingParserHelium" + lab.ExcitingVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.ExcitingVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """\s*\+-----------------------------------+\+\s*
\s*\|\s*EXCITING\s(?<version>helium\s*\S*) started\s*\|\s*
(?:\s*\|\sversion hash id:\s*(?<hashId>\S+)\s*\|)?""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/exciting/parser/parser-exciting/parser_exciting_helium.py",
    "--uri", "${mainFileUri}", "${mainFilePath}"),
  mainFileMatchPriority = 5,
  resList = Seq(
    "parser-exciting/parser_exciting_helium.py",
    "parser-exciting/exciting_parser_dos.py",
    "parser-exciting/exciting_parser_bandstructure.py",
    "parser-exciting/exciting_parser_input.py",
    "parser-exciting/setup_paths.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/exciting.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-exciting" -> "parsers/exciting/parser/parser-exciting/",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info",
    "python" -> "python-common/common/python/nomadcore"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
)