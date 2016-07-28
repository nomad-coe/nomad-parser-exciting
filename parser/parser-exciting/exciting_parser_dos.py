import xml.sax
import logging

class DosHandler(xml.sax.handler.ContentHandler):
    def __init__(self, backend):
        self.backend = backend
        self.dosSectionGIndex = -1
        self.inDos = False

    def endDocument(self):
        pass

    def startElement(self, name, attrs):
        if name == "dos":
            self.dosSectionGIndex = self.backend.openSection("x_exciting_section_dos")
            self.inDos = True
        elif name == "point" and self.inDos:
            self.backend.addValue("x_exciting_dos_value",float(attrs.getValue('dos')))
            self.backend.addValue("x_exciting_dos_energy",float(attrs.getValue('e')))
        # attrDict={}
        # for name in attrs.getNames():
        #     attrDict[name] = attrs.getValue(name)
        # logging.error("start element %s attr %s", name, attrDict)

    def endElement(self, name):
        if name == 'dos':
            self.inDos = False
            self.backend.closeSection("x_exciting_section_dos",self.dosSectionGIndex)
            self.dosSectionGIndex = -1
        # logging.error("end element %s", name)

    def startElementNS(self, name, qname, attrs):
        attrDict={}
        for name in attrs.getNames():
            attrDict[name] = attrs.getValue(name)
        logging.error("start element %s ns %s attr %s", name, qname, attrDict)

    def endElementNS(self, name, qname):
        logging.error("end element %s ns %s", name, qname)

    def characters(self, content):
        pass

def parseDos(inF, backend):
    handler = DosHandler(backend)
    logging.error("will parse")
    xml.sax.parse(inF, handler)
    logging.error("did parse")
