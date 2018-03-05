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
    logging.error("will parse")
    xml.sax.parse(inF, handler)
    logging.error("did parse")
