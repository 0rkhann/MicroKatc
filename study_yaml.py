"""Strict YAML loading for MicroKatc study files.

PyYAML follows YAML 1.1, which reads 1e-10 as text, NO/ON/yes/off as booleans, and keeps the
last of two duplicate keys silently. Study files need the opposite on all three counts.
"""

import os
import re

import yaml


class StudyError(Exception):
    """One or more problems in a study file; messages lists each with its location"""

    def __init__(self, messages):
        self.messages = list(messages)
        super().__init__("\n".join(self.messages))


class StrictLoader(yaml.SafeLoader):
    """SafeLoader with YAML 1.2-style floats, lowercase-only booleans and duplicate-key errors"""

    def construct_mapping(self, node, deep=False):
        lines = {}
        for key_node, _ in node.value:
            key = self.construct_object(key_node, deep=deep)
            line = key_node.start_mark.line + 1
            if key in lines:
                raise yaml.constructor.ConstructorError(
                    None,
                    None,
                    f"{key} appears twice (lines {lines[key]} and {line})",
                    key_node.start_mark,
                )
            lines[key] = line
        return super().construct_mapping(node, deep=deep)


# Copy the resolver table so SafeLoader itself is not changed
StrictLoader.yaml_implicit_resolvers = {
    first: [
        (tag, regexp)
        for tag, regexp in resolvers
        if tag not in ("tag:yaml.org,2002:bool", "tag:yaml.org,2002:float")
    ]
    for first, resolvers in yaml.SafeLoader.yaml_implicit_resolvers.items()
}
StrictLoader.add_implicit_resolver(
    "tag:yaml.org,2002:bool", re.compile(r"^(?:true|false)$"), list("tf")
)
StrictLoader.add_implicit_resolver(
    "tag:yaml.org,2002:float",
    re.compile(
        r"""^[-+]?(?:(?:[0-9][0-9_]*\.[0-9_]*|\.[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)?
        |[0-9][0-9_]*[eE][-+]?[0-9]+
        |\.(?:inf|Inf|INF)
        |\.(?:nan|NaN|NAN))$""",
        re.VERBOSE,
    ),
    list("-+0123456789."),
)


def load_yaml(path):
    """Reads a study file; raises StudyError with the file name and line for any YAML problem"""
    name = os.fspath(path)
    try:
        with open(name, encoding="utf-8") as f:
            doc = yaml.load(f, Loader=StrictLoader)
    except yaml.YAMLError as error:
        mark = getattr(error, "problem_mark", None) or getattr(
            error, "context_mark", None
        )
        where = f"{name} line {mark.line + 1}" if mark else name
        problem = getattr(error, "problem", None) or str(error)
        raise StudyError([f"{where}: {problem}"]) from error
    if not isinstance(doc, dict):
        raise StudyError(
            [
                f"{name}: the study must be a mapping of keys such as species, steps and conditions"
            ]
        )
    return doc
