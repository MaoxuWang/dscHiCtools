import setuptools

print (setuptools.find_packages())

setuptools.setup(
    name="dscHiCtools",
    version="1.1",
    author="Maoxu Wang",
    author_email="maoxuwang@stu.pku.edu.cn",
    description="A python package to process demultiplexed single cell Hi-C data",
    url="https://github.com/MaoxuWang/dscHiCtools/",
    packages=setuptools.find_packages(),
    entry_points={"console_scripts": ["dscHiCtools = dscHiCtools.__main__:main"]},
    include_package_data=True,

    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    install_requires=[
        "python-levenshtein>=0.12.0",
        "pysam>=0.16.0.1",
        "multiprocess>=0.70.6.1",
        "pytest==6.2.5",
        "pandas>=0.23.4",
        "pybktree==1.1",
    ],
    python_requires=">=3.6",
)
