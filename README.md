# Negative chi2 values from fit?

To build the framework and run it, follow the steps below:

```
cmsrel CMSSW_10_6_37
cd CMSSW_10_6_37/src
cmsenv

git clone https://github.com/P-H-Wagner/questions.git
cd questions/negChi2
scram b
cmsRun test/run_nano_rds_cfg.py
```
