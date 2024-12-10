import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="aip", # Replace with your own username
    version="0.0.0",
    author="M C Storer",
    author_email="mcs92@cam.ac.uk",
    description="A small example package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    entry_points={ 'console_scripts': [
        'filereader = filereader.__main__:main', 
        'filewriter = filewriter.__main__:main', 
        'aipAtomTypes = aipAtomTypes.__main__:main',
        'cmlgenerator = cmlgenerator.__main__:main', 
        'nwchemcmlutils = nwchemcmlutils.__main__:main']},
    include_package_data=True,
    zip_safe=False,
    data_files=[
          ('nwchemcmlutils/resources/',
           ['nwchemcmlutils/resources/camhpcparameters.json',
            'nwchemcmlutils/resources/ziggyparameters.json'])],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    #python_requires='>=3.6',
    python_requires='>=2.7',
)
