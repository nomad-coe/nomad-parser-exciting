/*
   Copyright 2016-2017 The NOMAD Developers Group

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 */
package eu.nomad_lab.parsers

import eu.{ nomad_lab => lab }
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import scala.collection.breakOut

object ExcitingParser extends SimpleExternalParserGenerator(
  name = "ExcitingParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("ExcitingParser")) ::
      ("parserId" -> jn.JString("ExcitingParser" + lab.ExcitingVersionInfo.version)) ::
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
  mainFileRe = """\s*=================================================+\s*
\s*\|\s*EXCITING\s(?<version>\S*) started\s*=
(?:\s*\|\sversion hash id:\s*(?<hashId>\S*)\s*=)?""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/exciting/parser/parser-exciting/parser_exciting.py",
    "--uri", "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-exciting/parser_exciting.py",
    "parser-exciting/exciting_parser_dos.py",
    "parser-exciting/exciting_parser_bandstructure.py",
    "parser-exciting/exciting_parser_gw.py",
    "parser-exciting/exciting_parser_input.py",
    "parser-exciting/exciting_parser_GS_input.py",
    "parser-exciting/exciting_parser_XS_input.py",
    "parser-exciting/exciting_parser_xs.py",
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
