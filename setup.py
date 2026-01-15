from setuptools import setup, find_packages

# with open("README.md", "r", encoding="utf-8") as fh:
#     long_description = fh.read()

setup(
    name="organelle_measure",
    version="0.0.1",
    author="Shixing Wang",
    author_email="wangshixing@wustl.edu",
    description="a toolbox to measure organlles",
    # long_description=long_description,
    # long_description_content_type="text/markdown",
    # url="https://github.com/pypa/sampleproject",
    # project_urls={
    #     "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    # },
    # classifiers=[
    #     "Programming Language :: Python :: 3",
    #     "License :: OSI Approved :: MIT License",
    #     "Operating System :: OS Independent",
    # ],
    # package_dir={"": "src"},
    # packages=setuptools.find_packages(where="src"),
    packages=find_packages(),
    python_requires=">=3.6",
)