package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object ExcitingParserSpec extends Specification {
  "ExcitingParserTest" >> {
    "test with json-events" >> {
      //No test file present at the moment; Replace the README when test file is present
      ParserRun.parse(ExcitingParser, "parsers/castep/test/examples/README.md", "json-events") must_== ParseResult.ParseSuccess
    }
  }

  "test with json" >> {
    ParserRun.parse(ExcitingParser, "parsers/exciting/test/examples/README.md", "json") must_== ParseResult.ParseSuccess
  }
}
