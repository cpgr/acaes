#include "acaesApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
acaesApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

acaesApp::acaesApp(InputParameters parameters) : MooseApp(parameters)
{
  acaesApp::registerAll(_factory, _action_factory, _syntax);
}

acaesApp::~acaesApp() {}

void 
acaesApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<acaesApp>(f, af, s);
  Registry::registerObjectsTo(f, {"acaesApp"});
  Registry::registerActionsTo(af, {"acaesApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
acaesApp::registerApps()
{
  registerApp(acaesApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
acaesApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  acaesApp::registerAll(f, af, s);
}
extern "C" void
acaesApp__registerApps()
{
  acaesApp::registerApps();
}
