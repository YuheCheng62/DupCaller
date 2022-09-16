import os
from setuptools import setup, find_packages


def read(filename):
    filepath = os.path.join(os.path.dirname(__file__), filename)
    file = open(filepath, 'r')
    return file.read()


setup(
    name='DupCaller',
    version='0.0.1',
    author='Yuhe Cheng',
    author_email='cheng.yuhe62@gmail.com',
    description='A somatic mutation caller for duplex sequencing data',
    long_description=read('README.md'),
    long_description_content_type="text/markdown",
    license='GNU General Public License v3.0',
    url='https://github.com/AlexandrovLab/DupCaller',
    package_dir = {'':'src'},
    packages=find_packages('src'),
    install_requires = [
        "numpy",
        "pysam",
        "biopython",
        "pandas"],
    scripts=["src/DupCallerTrim.py","src/DupCallerCall.py"],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Development Status :: Vanilla',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3'
    ],
)