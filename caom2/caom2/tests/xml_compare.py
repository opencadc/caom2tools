"""
NOTE: This is the xml_compare v1.0.5 module available on pypi
(https://pypi.python.org/pypi/xml_compare/1.0.5)
with a tiny patch on line 69 to make it work with Python 3.
The module hasn't been updated in almost 10 years,
hence the need to copy it over here and patch it.

Also patched the xml_compare function to deal with children nodes that are in
different order.

"""

import os

"""
This module implements a XML comparer

>>> from lxml import etree
>>> xml1 = "<a><b/><c>text</c><d/></a>"
>>> doc1 = etree.fromstring(xml1)
>>> xml2 = "<a><b/><c>text</c><d/></a>"
>>> doc2 = etree.fromstring(xml2)
>>> xml_compare(doc1,doc2)
True
>>> xml3 = "<a><b/><c>texto</c><d/></a>"
>>> doc3 = etree.fromstring(xml3)
>>> xml_compare(doc1,doc3)
False
>>> xml4 = "<a><b/><c> text </c><d/></a>"
>>> doc4 = etree.fromstring(xml4)
>>> xml_compare(doc1,doc4,reporter=sys.stderr.write,strip_whitespaces=True)
True
>>> xml5a = "<a>4567.23000</a>"
>>> xml5b = "<a>4567.23</a>"
>>> doc5a = etree.fromstring(xml5a)
>>> doc5b = etree.fromstring(xml5b)
>>> xml_compare(doc5a,doc5b,reporter=sys.stderr.write,strip_whitespaces=True,
                float_compare=True)
True
>>> xml6a = "<a v='4567.23'/>"
>>> xml6b = "<a v='4567.230000'/>"
>>> doc6a = etree.fromstring(xml6a)
>>> doc6b = etree.fromstring(xml6b)
>>> xml_compare(doc6a,doc6b,reporter=sys.stderr.write,strip_whitespaces=True,
                float_compare=True)
True
"""


def getNodePath(node):
    return node.getroottree().getpath(node)


def doStripWhitespaces(text):
    if text is None:
        return None
    else:
        return text.strip().replace('\n', '').replace('\r', '')


def text_compare(t1, t2, strip_whitespaces=False, float_compare=False):
    """ Text comparer for XML Text Nodes
    >>> text_compare("hi","hi")
    True
    >>> text_compare("hi","hi ")
    False
    >>> text_compare("hi","hi ", strip_whitespaces=True)
    True
    >>> text_compare("12.234","12.324", strip_whitespaces=True,
                     float_compare=True)
    False
    >>> text_compare("12.234","12.2340000", strip_whitespaces=True,
                     float_compare=True)
    True
    >>> text_compare("Hola","Algo", strip_whitespaces=True, float_compare=True)
    False
    """
    if not t1 and not t2:
        return True
    if float_compare:
        try:
            f1 = float(t1)
            f2 = float(t2)
            return f1 == f2
        except (ValueError, TypeError):
            pass
    if strip_whitespaces:
        return (doStripWhitespaces(t1) or '') == (doStripWhitespaces(t2) or '')
    else:
        return (t1 or '') == (t2 or '')


def doReport(reporter, x1, x2, errorMsg):
    if reporter:
        reporter(getNodePath(x1) + " " + getNodePath(
            x2) + os.linesep + errorMsg + os.linesep)


def xml_compare(x1, x2, reporter=None, strip_whitespaces=False,
                ignore_order=False, float_compare=False):
    if x1.tag != x2.tag:
        doReport(reporter, x1, x2,
                 'Tags do not match: %s and %s' % (x1.tag, x2.tag))
        return False

    for name, value in x1.attrib.items():
        if not text_compare(value, x2.attrib.get(name),
                            strip_whitespaces=strip_whitespaces,
                            float_compare=float_compare):
            doReport(reporter, x1, x2, 'Attributes do not match: %s=%r, %s=%r'
                     % (name, value, name, x2.attrib.get(name)))
            return False
    for name in x2.attrib.keys():
        if name not in x1.attrib:
            doReport(reporter, x1, x2,
                     'x2 has an attribute x1 is missing: %s' % name)
            return False
    if not text_compare(x1.text, x2.text, strip_whitespaces=strip_whitespaces,
                        float_compare=float_compare):
        doReport(reporter, x1, x2, 'text: %r != %r' % (x1.text, x2.text))
        return False

    if not text_compare(x1.tail, x2.tail, strip_whitespaces=strip_whitespaces,
                        float_compare=float_compare):
        doReport(reporter, x1, x2, 'tail: %r != %r' % (x1.tail, x2.tail))
        return False

    cl1 = x1.getchildren()
    cl2 = x2.getchildren()
    if len(cl1) != len(cl2):
        doReport(reporter, x1, x2,
                 'children length differs, %i != %i' % (len(cl1), len(cl2)))
        return False
    found = False
    for c1 in cl1:
        for c2 in cl2:
            if xml_compare(c1, c2, reporter,
                           strip_whitespaces=strip_whitespaces,
                           float_compare=float_compare):
                found = True
                break
        if not found:
            doReport(reporter, c1, c1,
                     'child %s not found in destination' % (c1.tag))
            return False
    return True


def _test():
    import doctest
    doctest.testmod()


if __name__ == "__main__":
    _test()
