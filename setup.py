from setuptools import setup

setup(
    name='pyEliasFano',
    version='0.0.1',
    description='pyEliasFano offers a **quasi-succinct** represenation for a monotone non-decreasing sequence of n integers from the universe [0 . . . m) occupying 2*n+n*ceil(log2(m/n)) bits.',
    url='https://github.com/rmrschub/pyEliasFano',
    author='René Schubotz',
    author_email='rene.schubotz@dfki.de',
    license='CC BY-NC-SA 4.0',
    packages=['pyEliasFano'],
    install_requires=['math',
                      'typing',
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