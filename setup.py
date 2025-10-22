from setuptools import setup, find_packages

setup(
        name="livermeth",
        version='0.1',
        packages=find_packages(),
        author='Manish Yadav',
        author_email='manish20072013 at gmail dot com',
        description='A package for analyzing liver methylation data',
        long_description=open('README.md', encoding='utf-8').read(),
        long_description_content_type='text/markdown',
        license='GPL3', 
        python_requires=">=3.10",
        
    )