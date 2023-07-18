from label_transfer_predictions import LabelTransferPredictor

predictor = LabelTransferPredictor(verbose=False)

predictor.predict_all(threshold=0.6)

predictor.plot_predictions_umap("Ruf_Zamojski_2021M")