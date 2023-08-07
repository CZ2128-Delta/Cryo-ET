from setuptools import find_packages, setup

setup(
    name="imodfile_parse_module",
    version="0.0",
    install_requires=[],
    packages=find_packages(),
    extras_require={
        "all": ["matplotlib", "pycocotools", "opencv-python", "onnx", "onnxruntime"],
        "dev": ["flake8", "isort", "black", "mypy"],
    },
)
