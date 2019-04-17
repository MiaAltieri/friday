import torch
from modules.python.models.inception import Inception3


class ModelHandler:
    @staticmethod
    def save_checkpoint(state, filename):
        torch.save(state, filename)

    @staticmethod
    def get_new_model(input_channels, num_classes=6):
        # get a new model
        model = Inception3(input_channels=input_channels, num_classes=num_classes,
                           aux_logits=False, transform_input=False)
        return model

    @staticmethod
    def load_optimizer(optimizer, checkpoint_path, gpu_mode):
        if gpu_mode:
            checkpoint = torch.load(checkpoint_path)
            optimizer.load_state_dict(checkpoint['encoder_optimizer'])
            for state in optimizer.state.values():
                for k, v in state.items():
                    if isinstance(v, torch.Tensor):
                        state[k] = v.cuda()
        else:
            checkpoint = torch.load(checkpoint_path, map_location='cpu')
            optimizer.load_state_dict(checkpoint['optimizer'])

        return optimizer

    @staticmethod
    def load_model_for_training(model_path, input_channels, num_classes):
        checkpoint = torch.load(model_path, map_location='cpu')
        epochs = checkpoint['epochs']
        model = ModelHandler.get_new_model(input_channels=input_channels, num_classes=num_classes)
        state_dict = checkpoint['state_dict']

        from collections import OrderedDict
        new_encoder_state_dict = OrderedDict()
        new_decoder_state_dict = OrderedDict()

        for k, v in state_dict.items():
            name = k
            if k[0:7] == 'module.':
                name = k[7:]  # remove `module.`
            new_encoder_state_dict[name] = v

        for k, v in state_dict.items():
            name = k
            if k[0:7] == 'module.':
                name = k[7:]  # remove `module.`
            new_decoder_state_dict[name] = v

        model.load_state_dict(new_encoder_state_dict)
        model.cpu()

        return model, epochs

