from xml.etree.ElementTree import XMLTreeBuilder
import sys


##Adi:  I was just playing around with this core code to see what are the ideas of this.
## you can delete this file with any problems

class orthoxmlParser(object):

    def __init__(self, outputFile):
        self.group_name = ''
        return

    def start(self, tag, attrib):
        print("start")
        print(tag)
        print(attrib)

    def end(self, tag):
        # Ignore closing tags
        print("end")
        pass
    def data(self, data):
        # Ignore data inside nodes
        pass
    def close(self):
        # Nothing special to do here
        return


target = orthoxmlParser(sys.stdout)
parser = XMLTreeBuilder(target=target)
with open('./test/simpleEx.orthoxml', 'rt') as f:
    for line in f:
        parser.feed(line)