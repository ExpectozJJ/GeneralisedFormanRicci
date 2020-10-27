from GeneralisedFormanRicci.frc import GeneralisedFormanRicci

def test_sample():
  data = [[0.8, 2.6], [0.2, 1.0], [0.9, 0.5], [2.7, 1.8], [1.7, 0.5], [2.5, 2.5], [2.4, 1.0], [0.6, 0.9], [0.4, 2.2]]
  for f in [0, 0.5, 1, 2, 3]:
      sc = GeneralisedFormanRicci(data, method = "rips", epsilon = f)
      sc.compute_forman()
