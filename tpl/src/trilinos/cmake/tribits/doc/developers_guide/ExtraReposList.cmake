TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES(
  ExtraRepo1  ""  GIT  someurl.com:/ExtraRepo1   ""          Continuous
  ExtraRepo2  packages/SomePackage/Blah  GIT  someurl2.com:/ExtraRepo2
                                                 NOPACKAGES  Nightly
  ExtraRepo3  ""  HG  someurl3.com:/ExtraRepo3   ""          Continuous
  ExtraRepo4  ""  SVN  someurl4.com:/ExtraRepo4  ""          Nightly
  )
