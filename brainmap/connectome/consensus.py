"""
conensus.py


Module for handling consensus data


"""

from lxml import etree

def convert_synapse_to_xml(data):
    root = etree.Element('conserved_synapses')
    for cell in data:
        xcell = etree.SubElement(root,'cell')
        xcell.set('name',cell)
        for cont in data[cell]:
            xcont = etree.SubElement(xcell,'contin')
            xcont.set('id',cont)
            xnum = etree.SubElement(xcont,'num_sections')
            xnum.text = str(data[cell][cont]['num_sections'])
            xsect = etree.SubElement(xcont,'sections')
            for s in data[cell][cont]['sections']:
                xname = etree.SubElement(xsect,'name')
                xname.text = s
            xpartner = etree.SubElement(xcont,'partners')
            for p in data[cell][cont]['partners']:
                xname = etree.SubElement(xpartner,'name')
                xname.text = p
            xneighbor = etree.SubElement(xcont,'neighbors')
            for n in data[cell][cont]['neighbors']:
                xname = etree.SubElement(xneighbor,'name')
                xname.text = n
            xcell.append(xcont)
        root.append(xcell)

    tree = etree.ElementTree(root)
    return tree      

def convert_xml_to_synapse(fin):
    data = {}
    tree = etree.parse(fin)
    root = tree.getroot()
    for cell in root.findall('cell'):
        cname = cell.get('name')
        data[cname] = {}
        cont = cell.findall('contin')
        for c in cont:
            contid = c.get('id')
            num = c.find('num_sections').text
            xsect = c.find('sections')
            sections = [s.text for s in xsect.findall('name')]
            xpartner = c.find('partners')
            partners = [p.text for p in xpartner.findall('name')]
            xneighbor = c.find('neighbors')
            neighbors = [p.text for p in xneighbor.findall('name')]
            data[cname][contid] = {'num_sections':int(num),
                                    'sections':sections,
                                    'partners':partners,
                                    'neighbors':neighbors}
    return data
                                    

