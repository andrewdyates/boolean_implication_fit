for line in open("fake_data.tab"):
  print ",".join(line.strip('\n\r').split('\t'))
