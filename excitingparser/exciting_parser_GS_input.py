# Copyright 2016-2018 The NOMAD Developers Group
# Copyright 2017-2018 Lorenzo Pardini
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import xml.sax
import logging
import numpy as np
from nomadcore.unit_conversion.unit_conversion import convert_unit_function
from nomadcore.unit_conversion.unit_conversion import convert_unit
from nomadcore.unit_conversion import unit_conversion


class InputHandler(xml.sax.handler.ContentHandler):
    def __init__(self, backend):
        self.backend = backend

    def startElement(self, name, attrs):
        if name == "libxc":        #libXC
            correlation = attrs.getValue("correlation")[3:]
            exchange = attrs.getValue("exchange")[3:]
            xcName = [correlation, exchange]
            for xc in xcName:
                gi = self.backend.openSection("section_XC_functionals")
                self.backend.addValue("XC_functional_name", xc)
                self.backend.closeSection("section_XC_functionals", gi)

    def endElement(self, name):
        pass

def parseInput(inF, backend):
    handler = InputHandler(backend)
    logging.info("will parse")
    xml.sax.parse(inF, handler)
    logging.info("did parse")
