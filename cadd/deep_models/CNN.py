from tensorflow.keras.layers import Conv2D, Input, MaxPool2D, Activation, GlobalMaxPooling2D, Dense, Flatten
from tensorflow.keras.models import Model

class CnnModel:
    def __init__(self):
        self.input_shape =

    def build_model(self):

        x_input = Input(self.input_shape)

        # first convolution
        x = Conv2D(64, (self.input_shape[0], 6), strides=(1, 1), padding='valid')(x_input)
        x = Activation('relu')(x)
        x = MaxPool2D((1, 2))(x)

        # second convolution
        x = Conv2D(64, (1, 6), strides=(1, 1), padding='valid')(x)
        x = Activation('relu')(x)
        x = MaxPool2D((1, 2))(x)

        # Global pooling
        x = GlobalMaxPooling2D()(x)

        # Dense
        x = Flatten()(x)
        x = Dense(64, activation='relu')(x)
        x = Dense(1, activation='sigmoid')(x)

        model = Model(x_input, outputs=x)
        return model

    def get_cnn_fingerprints(self):
        pass