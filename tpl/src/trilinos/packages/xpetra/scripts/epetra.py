import sys
sys.path.append('lib/')

import os
from string import Template
from ConfigParser import SafeConfigParser
from XpetraLib import *
from XpetraLibConfig import *


def buildFuncLineEpetra( functionNode ):

#TODO: clean up
    tree = etree.parse(conf_XMLclass)
    root = tree.getroot() # root == <doxygen>
    classNode = root[0]   # classNode == <compounddef>

    fullClassName = classNode.xpath('compoundname')[0].text # Tpetra::Map
    baseClassName = fullClassName.lstrip('Tpetra::')        # Map
    className = 'Epetra'+baseClassName                      # EpetraMap
##

    # <name> = function name
    name = functionNode.xpath('name')[0].text
    if name == baseClassName: name = className
    if name == '~'+baseClassName: name = '~'+className

    # <type> = return type of the function
    type = functionNode.xpath('type')[0].xpath("string()")

    # <argsstring>
    argsstring = functionNode.xpath('argsstring')[0].text

    # hack for Vector:
    # - add missing 'typename'
    # - do not add MultiVector inherited methods
    if 'magnitudeType' in type: type = 'typename ' + type
    if functionNode.xpath('//compoundname')[0].text == 'Tpetra::Vector':
        if name in ['replaceGlobalValue','sumIntoGlobalValue','replaceLocalValue','sumIntoLocalValue'] and 'size_t vectorIndex' in argsstring: return ''

    # <param> -> get list of arg name as a string 'GIDList, nodeIDList, LIDList'
    # Simple version
    #        paramList = functionNode.xpath('param/declname/text()')
    #        paramStr  = ', '.join(param for param in paramList)

    # More complete version
    paramStr  = ''
    paramNodes = functionNode.xpath('param')
    #print name
    for paramNode in paramNodes:
        n = paramNode.xpath('declname')[0].xpath("string()")
        if paramNode.xpath('type')[0].xpath("string()") in conf_TypeWrapped:
            paramStr += "toEpetra(" + n + ")"
        else:
            paramStr += n

        paramStr += ", "

    paramStr = paramStr.rstrip(', ')

    # briefdescription
    briefdescription = functionNode.xpath("briefdescription")[0].xpath("string()")

    if len(type) > 0:
        declStr = type + " " + name + argsstring
    else:
        declStr = name + argsstring
    declStr = declStr.rstrip().replace('typename ','')

    if 'TPETRA_DEPRECATED' in type: return ''
    if "const =0" in argsstring: return '' #hack for CrsMatrix

    # hack for MultiVector
    if name == "scale" and "Teuchos::ArrayView< const Scalar > alpha" in argsstring: return ''
    if name == "scale" and "const Scalar &alpha, const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &A" in argsstring: return ''

    # hack for Vector
    if name == 'EpetraVector' and 'Map' in argsstring: declStr = 'explicit ' + declStr
    if className == "EpetraVector" and 'ArrayView' in argsstring: return ''

    # hack for CrsMatrix
    if name == "EpetraCrsMatrix" and "const RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node> > &graph" in argsstring: return ''

    if name in conf_RemoveRefFunctionList: declStr = declStr.replace('&', '')

    descStr = "    //! " + briefdescription.lstrip().rstrip() + "\n"
    defStr  = "    " + declStr

    # implemented in .cpp
    if name in conf_FunctionInCppFile:
        defStr += ';'
        return descStr + defStr + "\n" + "\n";

    if name != className and name != "~"+className:
        defStr += " { "
        defStr += "XPETRA_MONITOR(\"" + className + "::" + name + "\"); "
        if len(type) > 0 and type != 'void': defStr += 'return '
        if type in conf_TypeWrapped: defStr += "toXpetra("

        if parser.has_option('replace',name):
            name = parser.get('replace', name)

        defStr += conf_memberName + "->" + name
        defStr += "(" + paramStr
        if type in conf_TypeWrapped: defStr += ")"
        defStr += "); }"

    # constructor
    if name == className:
        defStr += "\n      " + ": " + conf_memberName + "(Teuchos::rcp(new " + fullClassName.replace('Tpetra::','Epetra_')
        defStr += "(" + paramStr + "))) { }"

    # destructor
    if name == '~'+className:
        defStr += " { }"

    return descStr + defStr + "\n" + "\n";

####

xml_dir = trilinosRoot_dir + '/packages/tpetra/doc/xml/'
conf_dir = 'epetra/conf/'
tmpl_dir = 'epetra/tmpl/'
out_dir = '../src/'

for file in os.listdir(conf_dir):
    basename, extension = os.path.splitext(file)
    if extension == ".conf":
#### READ CONFIG ####
        parser = SafeConfigParser()
        parser.read(conf_dir + file)

        conf_XMLheaders = xml_dir + parser.get('io', 'XMLheaders')
        conf_XMLclass   = xml_dir + parser.get('io', 'XMLclass')
        conf_template   = tmpl_dir + parser.get('io', 'template')
        conf_output     = parser.get('io', 'output')

        conf_SkipFunctionList = set(parser.get('function', 'skip').split(';'))
        conf_RemoveRefFunctionList = set(parser.get('function', 'removeref').split(';'))
        if parser.has_option('function', 'inCppFile'):
            conf_FunctionInCppFile = set(parser.get('function', 'inCppFile').split(';'))
        else:
            conf_FunctionInCppFile = []
        conf_SkipHeaderList = set(parser.get('header', 'skip').split(';'))
        conf_memberName = parser.get('member', 'name')
        conf_TypeWrapped = set(parser.get('type', 'wrapped').split(';'))
#

        template = open(conf_template, 'r').read()
        out = Template(template)

        className = buildClassDefinition(conf_XMLclass, 'Epetra')
#unused        templateParam = buildTemplateParam2(conf_XMLclass)

        out = out.substitute(
            TMPL_HEADERS=buildHeader(className, 'epetra.py'),
            TMPL_INCLUDES=buildInclude(conf_XMLheaders, conf_SkipHeaderList),
            TMPL_TEMPLATE_PARAM=buildTemplateParam(conf_XMLclass),
            TMPL_CLASS=className,
            TMPL_INHERITANCE='  ' + parser.get('inheritance', 'parent').rstrip(),
            TMPL_DESTRUCTOR=buildDestructor(className),
            TMPL_PUBLIC_FUNCTIONS=buildClassFunctions(conf_XMLclass, conf_SkipFunctionList, buildFuncLineEpetra),
            TMPL_FOOTERS=buildFooter(className, short=False)
            )
        f = open(out_dir + conf_output, 'w')
        f.write(out)
        f.close()

