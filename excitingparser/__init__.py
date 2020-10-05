from .metainfo import m_env

from nomad.parsing.parser import FairdiParser
from excitingparser.exciting_parser import ExcitingOutput


class ExcitingParser(FairdiParser):
    def __init__(self):
        super().__init__(
            name='parsers/exciting', code_name='exciting', code_homepage='http://exciting-code.org/',
            mainfile_name_re=r'^.*.OUT(\.[^/]*)?$', mainfile_contents_re=(r'EXCITING.*started'))

    def parse(self, filepath, archive, logger=None):
        self._metainfo_env = m_env

        parser = ExcitingOutput(filepath, archive, logger)

        parser.parse()
