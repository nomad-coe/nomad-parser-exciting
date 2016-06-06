import xml.sax
import logging

class BandHandler(xml.sax.handler.ContentHandler):
    def __init__(self, backend):
        self.backend = backend
        self.bandSectionGIndex = -1
        self.inBand = False

    def endDocument(self):
        pass

    def startElement(self, name, attrs):
        if name == "bandstructure":
            self.bandSectionGIndex = self.backend.openSection("exciting_section_bandstructure")
            self.inBand = True
        elif name == "point" and self.inBand:
            self.backend.addValue("exciting_band_value",float(attrs.getValue('eval')))
            self.backend.addValue("exciting_band_k",float(attrs.getValue('distance')))
        # attrDict={}
        # for name in attrs.getNames():
        #     attrDict[name] = attrs.getValue(name)
        # logging.error("start element %s attr %s", name, attrDict)

    def endElement(self, name):
        if name == 'bandstructure':
            self.inBand = False
            self.backend.closeSection("exciting_section_bandstructure",self.bandSectionGIndex)
            self.bandSectionGIndex = -1
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

def parseBand(inF, backend):
    print("pippolo")
    handler = BandHandler(backend)
    logging.error("will parse")
    xml.sax.parse(inF, handler)
    logging.error("did parse")
