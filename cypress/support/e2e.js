before(() => {
    cy.visit("/")
    cy.get('#splash-button', {timeout: 40000}).should('be.enabled')
    cy.get('#splash-button').click()
})