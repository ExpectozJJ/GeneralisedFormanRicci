import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="GeneralisedFormanRicci", # Replace with your own username
    version="0.0.3",
    author="Wee JunJie",
    author_email="expectozjj@gmail.com",
    description="A class to compute the Generalised Forman-Ricci curvature for a Simplicial Complex from a given point cloud data.",
    long_description=long_description, 
    long_description_content_type="text/markdown",
    url="https://github.com/ExpectozJJ/GeneralisedFormanRicci",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)