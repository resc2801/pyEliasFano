from setuptools import setup

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
    long_description_content_type = "text/x-rst"
except(IOError, ImportError):
    long_description = open('README.md').read()
    long_description_content_type = 'text/markdown'

setup(
    name='pyEliasFano',
    version='0.0.2',
    description='pyEliasFano offers a **quasi-succinct** represenation for a monotone non-decreasing sequence of n integers from the universe [0 . . . m) occupying 2*n+n*ceil(log2(m/n)) bits.',
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    url='https://github.com/rmrschub/pyEliasFano',
    author='René Schubotz',
    author_email='rene.schubotz@dfki.de',
    license='CC BY-NC-SA 4.0',
    packages=['pyEliasFano'],
    install_requires=['typing',
                      'bitarray',
                      'unary_coding',
                      ],

    classifiers=[
#       'Development Status:: 3 - Alpha',
#        'License:: CC BY-NC-SA 4.0',
#        'Programming Language:: Python:: 3:: Only',
#        'Operating System:: OS Independent',
#        'Topic:: System:: Archiving:: Compression',
    ],
)