from setuptools import setup

setup(
    name='Spectator',
    author='Molly Peeples',
    author_email='molly@stsci.edu',
    url="https://github.com/peeples/spectator",
    version='0.1',
    py_modules=['drive_quick_look',
                'quick_look',
                'scrape_headers'
                'spectator'],  # Python modules to install (without the .py in the filename)
    scripts=['spectator.py']  # This is the full name of the script "beetle"; this will be installed to a bin/ directory
)
