
In order to make the fitpack module work, perform the following commands

   - cd stellapy.calculations/utils
   - rm *.so
   - f2py -c fitpack.pyf fitpack.f
   - cp fitpack.cpython-36m-x86_64-linux-gnu.so fitpack.so

Ofcourse in the last line <fitpack.cpython-36m-x86_64-linux-gnu.so> might be different on your computer, just make sure that the file created by f2py is renamed to fitpack.so.
