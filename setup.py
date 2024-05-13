"""
setup script
"""
import setuptools

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name="covizu",
    version="0.3",
    description="Visualization of SARS-CoV-2 diversity",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Art Poon, Emmanuel Wong, Kaitlyn Wade, Molly Liu, Roux-Cil Ferreira, "
           "Laura Munoz Baena, Abayomi Olabode",
    url="https://github.com/PoonLab/covizu",
    packages=setuptools.find_packages(),
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_data={'covizu': [
        'data/NC_045512.fa',
        'data/ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf',
        'data/pango-designation/lineages.csv'
    ]}
)
