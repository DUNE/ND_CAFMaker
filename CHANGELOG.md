# Changelog
[[Format loosely based on <https://keepachangelog.com/en/0.3.0>]]

##### current

##### [v4.7.0] -- 2024-10-16
* Add TMS Branch reco with updated StandardRecord format with one entry per spill and trigger logic ([PR #87](https://github.com/DUNE/ND_CAFMaker/pull/87))
* Add track matching logic for 2x2 reconstructions (Pandora, SPINE to MINERvA) ([PR #86](https://github.com/DUNE/ND_CAFMaker/pull/86))
* Add trigger time in meta information of SPINE and MINERvA ([PR #84](https://github.com/DUNE/ND_CAFMaker/pull/84))
* Fix: TMS start and end z position for tracks ([PR #83](https://github.com/DUNE/ND_CAFMaker/pull/83))

##### [v4.6.5] -- 2024-09-04
* Fix: Trigger matching second to nanosecond conversion ([PR #82](https://github.com/DUNE/ND_CAFMaker/pull/82))

##### [v4.6.4] -- 2024-09-03
* Changed ML-Reco branch filler to read new SPINE format ([PR #74](https://github.com/DUNE/ND_CAFMaker/pull/74))
* Added POT information for Data from file  ([PR #78](https://github.com/DUNE/ND_CAFMaker/pull/78))
* Added flash info in ML-Reco branch filler  ([PR #80](https://github.com/DUNE/ND_CAFMaker/pull/80))
* Added the possibility to read Truth info from Edepsim file  ([PR #77](https://github.com/DUNE/ND_CAFMaker/pull/77))

##### [v4.6.3] -- 2024-05-17
* Fix: Fill the flat CAFs correctly ([PR #65](https://github.com/DUNE/ND_CAFMaker/pull/65))
* Fix: Correction to the truth backtracking for ML-Reco branch filler ([PR #68](https://github.com/DUNE/ND_CAFMaker/pull/68))

##### [v4.6.2] -- 2024-04-05
* Fix: Added some failsafe when MINERvA Trajectory number overloads its capacity ([PR #59](https://github.com/DUNE/ND_CAFMaker/pull/59))
* Fix: Some correction to the truth backtracking ([PR #61](https://github.com/DUNE/ND_CAFMaker/pull/61))

##### [v4.6.1] -- 2024-02-13
* Fix: actually use GENIE true vertex info as intended ([PR #56](https://github.com/DUNE/ND_CAFMaker/pull/56))
* Fix: don't use invalidated iterators in trigger group loop ([PR #57](https://github.com/DUNE/ND_CAFMaker/pull/57))

##### [v4.6.0] -- 2024-01-29
* Stop using ML-reco passthrough truth to set locations ([PR #55](https://github.com/DUNE/ND_CAFMaker/pull/55))

##### [v4.5.0] -- 2024-01-20
* Placeholder for POT setting ([PR #52](https://github.com/DUNE/ND_CAFMaker/pull/52))

##### [v4.4.1] -- 2023-12-21
* Fix truth filling for low-level reco ([PR #49](https://github.com/DUNE/ND_CAFMaker/pull/49))

##### [v4.4.0] -- 2023-11-19
* Bump to duneanaobj version v3_02_01 ([PR #46](https://github.com/DUNE/ND_CAFMaker/pull/46))
* Add licensing info ([PR #45](https://github.com/DUNE/ND_CAFMaker/pull/45))
* Full MINERvA branch filling (([PR #39](https://github.com/DUNE/ND_CAFMaker/pull/39)))
* Match triggers between MINERvA & ML reco ([PR #42](https://github.com/DUNE/ND_CAFMaker/pull/42))
* Support running with pass-through truth only ([PR #41](https://github.com/DUNE/ND_CAFMaker/pull/41))

##### [v4.3.0] -- 2023-11-09
* Use true trajectory ID passed through ML-reco to avoid some hacks [PR #38](https://github.com/DUNE/ND_CAFMaker/pull/38)
* Begin propagating trigger information through ML reco [PR #34](https://github.com/DUNE/ND_CAFMaker/pull/34)
* Fix lepton 4-momentum filling [PR #31](https://github.com/DUNE/ND_CAFMaker/pull/31)

##### [v4.2.0] -- 2023-10-30
* Use fully propagated `vertexID` information from upstream to do correct truth matching
  ([PR #28](https://github.com/DUNE/ND_CAFMaker/pull/28))
* Minor bugfix to file opening ([PR #32](https://github.com/DUNE/ND_CAFMaker/pull/32))

##### [v4.1.0] -- 2023-10-03
* ML-Reco filling updates ([PR #26](https://github.com/DUNE/ND_CAFMaker/pull/26), [PR #29](https://github.com/DUNE/ND_CAFMaker/pull/29))
* Truth-matching updates ([PR #27](https://github.com/DUNE/ND_CAFMaker/pull/27))

##### [v4.0.0] -- 2023-08-21
* Major restructuring ([PR #25](https://github.com/DUNE/ND_CAFMaker/pull/25))
  * Emit one `StandardRecord` per _reco trigger_, rather than per GENIE event
  * Introduce `TruthMatcher` to match truth info between upstream GENIE and truth passed through reco
  * Major rewrite of ND-LAr ML-reco filling

##### [v3.6.0] -- 2023-07-17
* Bring back flat-CAFs and actually enable them ([PR #13](https://github.com/DUNE/ND_CAFMaker/pull/13))

##### [v3.5.0] -- 2023-05-23  ([PR #23](https://github.com/DUNE/ND_CAFMaker/pull/23))
* Rewrite `NDLArMLRecoBranchFiller` to read "ana-lite" files output by official ML Reco tools 

##### [v3.4.2] -- 2023-03-28
* Improve documentation, minor cleanups

##### [v3.4.0] -- 2022-08-25
* Updates to TMS ([PR #9](https://github.com/DUNE/ND_CAFMaker/pull/9))
* Fill in SRTrack::dir for ND-LAr tracks ([PR #15](https://github.com/DUNE/ND_CAFMaker/pull/15))
* Update to e20 build ([PR #16](https://github.com/DUNE/ND_CAFMaker/pull/16))
* Eliminate "dumptree" ([PR #12](https://github.com/DUNE/ND_CAFMaker/pull/12))

##### [v3.3.1] -- 2022-08-18
* Minor fixes

##### [v3.3.0] -- 2022-08-18
* Include SAND

##### [v3.2.0] -- 2023-07-04
* More production tooling  ([PR #7](https://github.com/DUNE/ND_CAFMaker/pull/7))

##### [v3.1.0] -- 2022-06-21
* Production tooling ([PR #6](https://github.com/DUNE/ND_CAFMaker/pull/6))

##### [v3.0.0] -- 2022-05-24
* First version with reco: introducing ND-LAr reco ([PR #5](https://github.com/DUNE/ND_CAFMaker/pull/5))

##### [v2.2.0] -- 2021-07-28

##### [v2.1.0] -- 2021-01-15

##### [v2.0.0] -- 2020-05-07

##### [v1.0.0] -- 2020-04-29