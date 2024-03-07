//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "acaesTestApp.h"
#include "acaesApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
acaesTestApp::validParams()
{
  InputParameters params = acaesApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

acaesTestApp::acaesTestApp(InputParameters parameters) : MooseApp(parameters)
{
  acaesTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

acaesTestApp::~acaesTestApp() {}

void
acaesTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  acaesApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"acaesTestApp"});
    Registry::registerActionsTo(af, {"acaesTestApp"});
  }
}

void
acaesTestApp::registerApps()
{
  registerApp(acaesApp);
  registerApp(acaesTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
acaesTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  acaesTestApp::registerAll(f, af, s);
}
extern "C" void
acaesTestApp__registerApps()
{
  acaesTestApp::registerApps();
}
