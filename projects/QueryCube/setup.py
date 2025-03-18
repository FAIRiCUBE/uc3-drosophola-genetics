from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.0.1'
DESCRIPTION = 'Basic Package for CodeSharing'
LONG_DESCRIPTION = 'A package that allows to .....'

# Setting up
setup(
    name="WormPickerOOP",
    version=VERSION,
    author="Sonja Steindl",
    author_email="<sonja.steindl@nhm-wien.ac.at>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    install_requires=[],
    keywords=['python', 'video', 'stream', 'video stream', 'camera stream', 'sockets'],
    classifiers=[
        "Development Status :: 1 ",
        "Intended Audience ::  ",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix"
    ]
)