import os

with open(r'../he3.h', 'r') as header:
    with open(r'he3lib.py', 'w') as lib:
        lib.write("""'''
he3lib.py

Python interface to he3lib, a library of liquid helium-3 properties
'''
import numpy, ctypes

he3lib = ctypes.cdll.LoadLibrary('%s/libhe3.so')
""" % os.path.split(os.getcwd())[0])
        for line in header:
            try:
                line = line.strip()
                try:
                    docstr = '"' + line.split('/*')[1].split('*/')[0].strip().replace('"', "''") + '"; '
                    line = line.split('/*')[0]
                except:
                    docstr = ''
                if line.startswith('extern double'):
                    const_name = line.split('extern double')[1].split(';')[0].strip()
                    python_name = const_name.rstrip('_')
                    lib.write('''
def %s(): %sreturn ctypes.c_double.in_dll(he3lib, "%s").value
%s_const = ctypes.c_double.in_dll(he3lib, "%s").value
''' % (python_name, docstr, const_name, python_name, const_name))
                elif line.startswith('double') and line.count('(') == 1 and line.count(')') == 1:
                    func_name = line.split('double')[1].split('(')[0].strip()
                    python_name = func_name.rstrip('_')
                    args = [arg.strip() for arg in
                        line.split('(')[1].split(')')[0].strip().replace('double *', '').split(',')]
                    for arg in args:
                        if len(arg) < 1 or not arg[0].isalpha() or not arg.isalnum():
                            raise ValueError('unsupported argument list')
                    lib.write('''
he3lib.%s.restype = ctypes.c_double
def %s_scalar(%s): %sreturn he3lib.%s(%s)
%s = numpy.frompyfunc(%s_scalar, %d, 1)
''' % (func_name, python_name, ', '.join(args), docstr, func_name,
    ', '.join(['ctypes.byref(ctypes.c_double(%s))' % arg for arg in args]),
    python_name, python_name, len(args)))
            except Exception, ex:
                print '%s: %s' % (line, ex)
