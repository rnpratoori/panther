#include "pantherApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
pantherApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

pantherApp::pantherApp(InputParameters parameters) : MooseApp(parameters)
{
  pantherApp::registerAll(_factory, _action_factory, _syntax);
}

pantherApp::~pantherApp() {}

void
pantherApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<pantherApp>(f, af, s);
  Registry::registerObjectsTo(f, {"pantherApp"});
  Registry::registerActionsTo(af, {"pantherApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
pantherApp::registerApps()
{
  registerApp(pantherApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
pantherApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  pantherApp::registerAll(f, af, s);
}
extern "C" void
pantherApp__registerApps()
{
  pantherApp::registerApps();
}
