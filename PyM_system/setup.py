from setuptools import setup

setup(name='PyM',
      version='2006',
      description='Py Mathematics (PyM)',
      url='https://mat-web.upc.edu/people/sebastia.xambo/PyECC.html',
      author='N. Sayols & S. Xambó',
      author_email='sebastia.xambo@upc.edu',
      #license='MIT',
      packages=['PyM','PyM.PyWIT','PyM.PyWIT.PyECC','PyM.PyWIT.PyECC.arithmetic'],
      install_requires=[
          'numpy'
      ],
      zip_safe=False)