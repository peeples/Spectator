from setuptools import setup

setup(
    name='Spectator',
    version='0.1',
    packages=['spectator'],
    url='https://github.com/peeples/Spectator',
    license='',
    author='Molly Peeples',
    author_email='molly@stsci.edu',
    description='',
    install_requires=[
        "numpy",
        "scipy",
        "astropy",
        "fitsio",
        "matplotlib",
    ],
    entry_points={
        'console_scripts': [
            'spectator = spectator.spectator:main'
        ]
    }
)
