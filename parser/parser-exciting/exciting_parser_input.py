import xml.sax
import logging

class InputHandler(xml.sax.handler.ContentHandler):
    def __init__(self, backend):
        self.backend = backend
        self.inputSectionGIndex = -1
        self.inInput = False

    def endDocument(self):
        pass

    def startElement(self, name, attrs):
        if name == "structure":
            self.inputSectionGIndex = self.backend.openSection("section_system")
            self.inInput = True
        elif name == "atom" and self.inInput:
            g = attrs.getValue('coord')
            print ("coord=", g)
#            self.backend.addValue("atom_positions",float(attrs.getValue('coord')))
#            self.backend.addValue("x_exciting_dos_energy",float(attrs.getValue('e')))
        # attrDict={}
        # for name in attrs.getNames():
        #     attrDict[name] = attrs.getValue(name)
        # logging.error("start element %s attr %s", name, attrDict)

    def endElement(self, name):
        if name == 'structure':
            self.inInput = False
            self.backend.closeSection("section_system",self.inputSectionGIndex)
            self.inputSectionGIndex = -1
        # logging.error("end element %s", name)

#    def startElementNS(self, name, qname, attrs):
#        attrDict={}
#        for name in attrs.getNames():
#            attrDict[name] = attrs.getValue(name)
#        logging.error("start element %s ns %s attr %s", name, qname, attrDict)

#    def endElementNS(self, name, qname):
#        logging.error("end element %s ns %s", name, qname)

#    def characters(self, content):
#        pass

def parseInput(inF, backend):
    handler = InputHandler(backend)
    logging.error("will parse")
    xml.sax.parse(inF, handler)
    logging.error("did parse")
