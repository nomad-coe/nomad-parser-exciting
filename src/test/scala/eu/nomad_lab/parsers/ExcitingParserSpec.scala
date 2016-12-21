package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object ExcitingParserSpec extends Specification {
  "ExcitingParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(ExcitingParser, "parsers/exciting/test/examples/TiO2/INFO.OUT", "json-events") must_== ParseResult.ParseSuccess
    }
  }

  "test with json" >> {
    ParserRun.parse(ExcitingParser, "parsers/exciting/test/examples/TiO2/INFO.OUT", "json") must_== ParseResult.ParseSuccess
  }

  "ExcitingParserHeliumTest" >> {
    "test with json-events" >> {
      ParserRun.parse(ExcitingParserHelium, "parsers/exciting/test/examples/helium/INFO.OUT", "json-events") must_== ParseResult.ParseSuccess
    }
  }

  "test with json" >> {
    ParserRun.parse(ExcitingParserHelium, "parsers/exciting/test/examples/helium/INFO.OUT", "json") must_== ParseResult.ParseSuccess
  }
}
