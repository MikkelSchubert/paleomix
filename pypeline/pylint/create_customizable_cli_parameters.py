import logilab.astng

from logilab.astng import MANAGER

def _disable_infer(_self, *_args, **kwargs):
    raise logilab.astng.InferenceError()


def hashlib_transform(module):
    for cls_obj in module.get_children():
        if not isinstance(cls_obj, logilab.astng.Class):
            continue

        for member in cls_obj.get_children():
            if not isinstance(member, logilab.astng.Function):
                continue

            try:
                if 'pypeline.atomiccmd.builder.create_customizable_cli_parameters' in member.decoratornames():
                    if member.type == 'method':
                        member.type = 'classmethod'


                    # Crude workaround to spurious errors:
                    # Prevent pylint from attempting to infer the return type
                    member.infer_call_result = _disable_infer

            except logilab.astng.UnresolvableName:
                pass


def register(*_args, **_kwargs):
    MANAGER.register_transformer(hashlib_transform)
