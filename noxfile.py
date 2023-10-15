import nox

SOURCES = (
    "noxfile.py",
    "paleomix",
    "setup.py",
    "tests",
)


nox.options.sessions = [
    "style",
    "tests",
]


@nox.session
def style(session: nox.Session) -> None:
    session.install("black~=23.9.1")
    session.run("black", "--check", *SOURCES)
    session.install("isort~=5.12.0")
    session.run("isort", "--check-only", *SOURCES)


@nox.session
def lints(session: nox.Session) -> None:
    session.install("ruff==0.0.292")
    session.run("ruff", *SOURCES)


@nox.session()
def typing(session: nox.Session) -> None:
    session.install(".")
    session.install("pytest~=7.4")
    session.install("nox~=2023.4.22")
    session.install("pyright==1.1.330")
    session.run("pyright", *SOURCES)


@nox.session()
def tests(session: nox.Session) -> None:
    session.install("-e", ".")
    session.install("pytest~=7.4")
    session.install("coverage[toml]~=7.3")
    session.install("pytest-cov~=4.1")

    # TODO: Treat warnings as errors
    # "-Werror",
    session.run(
        "python3",
        "-X",
        "dev",
        "-m",
        "pytest",
        "tests",
        "--cov",
        "paleomix",
        "--cov",
        "tests",
        "--cov-report=xml",
        "--cov-report=term-missing",
        "-m",
        "not slow",
        "--quiet",
    )


@nox.session(python=["3.7", "3.8", "3.9", "3.10", "3.11", "3.12"])
def full_tests(session: nox.Session) -> None:
    session.install(".")
    session.install("pytest~=7.4")

    # TODO: add "-Werror"
    session.run(
        "python3",
        "-X",
        "dev",
        "-m",
        "pytest",
        "tests",
        "--quiet",
    )
