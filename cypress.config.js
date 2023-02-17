module.exports = {
  // These settings apply everywhere unless overridden
  pageLoadTimeout : 10000,
  requestTimeout : 10000,
  defaultCommandTimeout: 10000,
  viewportWidth: 1920,
  viewportHeight: 1080,
  e2e: {
    testIsolation: false,
    baseUrl: 'http://localhost:8001',
    setupNodeEvents(on, config) {
      // implement node event listeners here
    },
  },
};
